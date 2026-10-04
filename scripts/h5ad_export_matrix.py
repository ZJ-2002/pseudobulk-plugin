# 中文注释：h5ad X 层导出（2026-10-04 节点审计怪癖 4c 的补位节点）。
# 输入构念：任意 h5ad——引擎自产（X 为 lzf 压缩 CSR 组，h5py 原生识别
# filter 32000）或宿主/上游 anndata 写的（gzip/CSR 或 dense）均可。
# backed='r' 按行块切片，内存有界；导出为长表 parquet（cell_id/gene/value
# 三列，行组=chunk 个细胞），供 SQL/DataFusion/parquet_sql 消费——引擎
# 原生节点只有 obs 层工厂，X 层此前无任何导出路径。
# 稀疏矩阵只导非零元（CSR 语义，0=未检测/未表达，与缺失在下游不可区分，
# 记入 log）；dense 矩阵同样只导非零元以防空问爆炸，语义同上。
# 值一律 cast float64（整数计数 ≤2^53 精确）。
# 失败处理：文件非 h5ad/层不存在/行切片不被支持且矩阵放不进内存 → 退出 2。
# 解释边界：纯忠实搬运，不归一化不过滤；X 的语义（原始计数 vs 标准化值）
# 由上游决定，log 记录 dtype/稀疏结构/非零元数供下游判断。
import os
import sys

import anndata as ad
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

LAYER = os.environ.get("EXPORT_LAYER", "X").strip()
CHUNK = int(os.environ.get("EXPORT_CHUNK", "2000"))
INPUT = os.environ.get("AUTONOMICS_INPUT0")
OUT_MATRIX = os.environ["AUTONOMICS_OUTPUT0"]
OUT_OBS = os.environ["AUTONOMICS_OUTPUT1"]
OUT_VAR = os.environ["AUTONOMICS_OUTPUT2"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT3"]

log_lines = []


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


try:
    adata = ad.read_h5ad(INPUT, backed="r")
except Exception as error:  # noqa: BLE001 — 打不开即失败，无重试语义
    fail(f"cannot open h5ad: {error}")

if LAYER and LAYER != "X":
    if LAYER not in adata.layers:
        fail(f"layer missing: {LAYER}")
    matrix = adata.layers[LAYER]
    log_lines.append(f"layer\t{LAYER}")
else:
    matrix = adata.X
    log_lines.append("layer\tX")

n_cells, n_genes = adata.n_obs, adata.n_vars
obs_names = adata.obs_names.astype(str)
var_names = adata.var_names.astype(str)
log_lines.append(f"shape\tcells={n_cells}\tgenes={n_genes}\tdtype={getattr(matrix, 'dtype', '?')}")

# 行切片探针：backed CSR 支持按行切块；不支持（dense backed/部分层）时
# 全量读入（超内存会自然失败并给出可读错误）。
try:
    probe = matrix[0:min(CHUNK, n_cells)]
    probe = None  # 立即释放；探针只为触发 backed 路径的异常
    row_slicing = True
except Exception:
    matrix = matrix[:] if hasattr(matrix, "__getitem__") else matrix
    row_slicing = False
    if sp.issparse(matrix):
        matrix = matrix.tocsr()
    else:
        matrix = np.asarray(matrix)
    log_lines.append("WARN row slicing unsupported in backed mode; matrix loaded fully")


def row_block(i, j):
    return matrix[i:j]


schema = pa.schema([
    ("cell_id", pa.string()),
    ("gene", pa.string()),
    ("value", pa.float64()),
])

n_entries = 0
writer = pq.ParquetWriter(OUT_MATRIX, schema, compression="zstd")
try:
    for start in range(0, n_cells, CHUNK):
        stop = min(start + CHUNK, n_cells)
        block = row_block(start, stop)
        if sp.issparse(block):
            block = block.tocsr()
            block.eliminate_zeros()
            coo = block.tocoo()
            rows = np.asarray(coo.row)
            cols = np.asarray(coo.col)
            values = np.asarray(coo.data, dtype=np.float64)
        else:
            dense = np.asarray(block, dtype=np.float64)
            rows_local, cols_local = np.nonzero(dense)
            rows = rows_local
            cols = cols_local
            values = dense[rows_local, cols_local]
        block_cells = obs_names[start:stop]
        table = pa.table({
            "cell_id": np.asarray(block_cells, dtype=object)[rows],
            "gene": np.asarray(var_names, dtype=object)[cols],
            "value": values,
        })
        writer.write_table(table)
        n_entries += len(values)
        print(f"rows {start}-{stop}/{n_cells} nnz_cum={n_entries}", file=sys.stderr)
finally:
    writer.close()

adata.obs.to_csv(OUT_OBS, sep="\t", index_label="cell_id")
adata.var.to_csv(OUT_VAR, sep="\t", index_label="gene_id")

log_lines.append(
    f"entries\t{n_entries}\t"
    "note: nonzero-only export; zero vs missing indistinguishable downstream"
)
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write("\n".join(log_lines) + "\n")

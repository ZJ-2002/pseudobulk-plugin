# 中文注释：供者×细胞类型 pseudobulk 聚合。
# 输入构念：细粒度 scRNA/snRNA h5ad；obs[sample_col]=供者/动物（生物学单位），
# obs[celltype_col]=注释。规则来源：方案 §15.3——统计单位是供者，不按细胞数
# 拆伪重复；同供者多库天然被求和合并（库水平设计须在上游单独处理）。
# 失败处理：列缺失/矩阵空 → 退出 2；非整数矩阵 → WARN + 原值求和，解释边界
# （可能是标准化值，不是原始计数）写入 log 与 counts 表头注释行之外的正文字段。
# 解释边界：本节点不做 QC、不做批次校正、不注释细胞类型；输出仅是计数层。
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SAMPLE_COL = os.environ["PB_SAMPLE_COL"]
CELLTYPE_COL = os.environ["PB_CELLTYPE_COL"]
MIN_CELLS = int(os.environ.get("PB_MIN_CELLS", "10"))
LAYER = os.environ.get("PB_LAYER", "").strip()
INPUT = os.environ.get("AUTONOMICS_INPUT0")
OUT_COUNTS = os.environ["AUTONOMICS_OUTPUT0"]
OUT_SAMPLES = os.environ["AUTONOMICS_OUTPUT1"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT2"]

log_lines = []


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


adata = ad.read_h5ad(INPUT)
log_lines.append(f"input\t{INPUT}\tcells={adata.n_obs}\tgenes={adata.n_vars}")
for column in (SAMPLE_COL, CELLTYPE_COL):
    if column not in adata.obs.columns:
        fail(f"obs column missing: {column}")

if LAYER:
    if LAYER not in adata.layers:
        fail(f"layer missing: {LAYER}")
    matrix = adata.layers[LAYER]
    log_lines.append(f"layer\t{LAYER}")
else:
    matrix = adata.X
    log_lines.append("layer\tX")

if sp.issparse(matrix):
    matrix = matrix.tocsr()
    dense_probe = matrix.data[:10000] if matrix.nnz else np.array([])
else:
    matrix = sp.csr_matrix(np.asarray(matrix))
    dense_probe = matrix.data
is_integer = dense_probe.size == 0 or np.allclose(dense_probe, np.round(dense_probe))
if not is_integer:
    log_lines.append(
        "WARN matrix values are not integers; aggregating as-is. "
        "If values were normalised downstream must NOT model them as raw counts (plan §15.3)."
    )

# 分组键：供者 × 细胞类型；转 str 防止 categorical 排序歧义。
samples = adata.obs[SAMPLE_COL].astype(str).to_numpy()
celltypes = adata.obs[CELLTYPE_COL].astype(str).to_numpy()
keys = pd.Series([f"{s}__{c}" for s, c in zip(samples, celltypes)])

sizes = keys.value_counts()
kept = sorted(sizes[sizes >= MIN_CELLS].index)
dropped = sizes[sizes < MIN_CELLS]
for key, count in dropped.sort_index().items():
    log_lines.append(f"dropped\t{key}\tn_cells={count}\t<min_cells={MIN_CELLS}")
if not kept:
    fail(f"no pseudobulk group reaches min_cells={MIN_CELLS}")

# one-hot 稀疏乘法求和：cells×genes 乘 (groups×cells)^T → groups×genes。
index = {key: i for i, key in enumerate(kept)}
rows = keys.map(index).to_numpy()
mask = ~pd.isna(rows)
row_idx = rows[mask].astype(int)
col_idx = np.arange(matrix.shape[0])[mask]
onehot = sp.csr_matrix(
    (np.ones(len(row_idx), dtype=np.float64), (row_idx, col_idx)),
    shape=(len(kept), matrix.shape[0]),
)
aggregated = onehot @ matrix.astype(np.float64)
aggregated = aggregated.toarray()

genes = np.asarray(adata.var_names, dtype=str)
counts = pd.DataFrame(aggregated.T, index=genes, columns=kept)
if is_integer:
    counts = counts.astype(np.int64)
counts.index.name = "gene"
counts.to_csv(OUT_COUNTS, sep="\t")

sample_rows = []
for key in kept:
    sample, celltype = key.split("__", 1)
    sample_rows.append((key, sample, celltype, int(sizes[key])))
pd.DataFrame(
    sample_rows, columns=["pseudobulk_id", "sample", "celltype", "n_cells"]
).to_csv(OUT_SAMPLES, sep="\t", index=False)

log_lines.append(f"groups_kept={len(kept)}\tgenes={len(genes)}\tinteger={is_integer}")
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write("\n".join(log_lines) + "\n")

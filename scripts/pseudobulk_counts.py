# 中文注释：供者×细胞类型 pseudobulk 聚合。
# 输入构念：细粒度 scRNA/snRNA h5ad；obs[sample_col]=供者/动物（生物学单位），
# obs[celltype_col]=注释。规则来源：方案 §15.3——统计单位是供者，不按细胞数
# 拆伪重复；同供者多库天然被求和合并（库水平设计须在上游单独处理）。
# 失败处理（2026-10-05 审计 F07 收口）：
#   列缺失/矩阵空/基因名重复/负值/非有限值 → 退出 2；
#   非整数矩阵 → 默认硬错误（allow_non_integer=true 才放行，WARN 明示进 log）；
#   整数性/非负/有限校验覆盖稀疏矩阵**全部**存储条目（不再只探前 1 万个）；
#   obs 供者/细胞类型缺失值 → 显式剔除进 log（缺失率 100% 才报错），不静默并入。
# 分组键：供者×细胞类型按 (sample, celltype) 元组分组（pandas MultiIndex 语义），
# 不再用 "sample__celltype" 字符串拼接做组键——拼接键会把 ("a__b","c") 与
# ("a","b__c") 合并成同一组。输出 ID 仍为拼接串；若不同元组拼出相同串，
# 改用百分号编码（"_"→%5F）保证 ID 单射，log 记 id_encoding。
# schema：默认 deseq2——counts 首列 gene_id、metadata 首列 sample_id 且含
# condition 列（condition_column 参数从 obs 取，组内必须恒定），表头直接对上
# DESeq2/limma-voom；legacy 保留旧表头（gene / pseudobulk_id，无 condition）。
# 资源边界：聚合保持稀疏（onehot @ matrix），导出按基因分块稠密化写 TSV，
# 峰值内存受 块基因数×组数 控制，与 细胞数×基因数 全量稠密化脱钩。
# 解释边界：本节点不做 QC、不做批次校正、不注释细胞类型；输出仅是计数层。
import os
import sys
from urllib.parse import quote

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SAMPLE_COL = os.environ["PB_SAMPLE_COL"]
CELLTYPE_COL = os.environ["PB_CELLTYPE_COL"]
MIN_CELLS = int(os.environ.get("PB_MIN_CELLS", "10"))
LAYER = os.environ.get("PB_LAYER", "").strip()
SCHEMA = os.environ.get("PB_SCHEMA", "deseq2").strip().lower()
ALLOW_NON_INTEGER = os.environ.get("PB_ALLOW_NON_INTEGER", "false").strip().lower() == "true"
CONDITION_COL = os.environ.get("PB_CONDITION_COLUMN", "").strip()
INPUT = os.environ.get("AUTONOMICS_INPUT0")
OUT_COUNTS = os.environ["AUTONOMICS_OUTPUT0"]
OUT_SAMPLES = os.environ["AUTONOMICS_OUTPUT1"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT2"]

log_lines = []


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


if SCHEMA not in ("deseq2", "legacy"):
    fail(f"unsupported schema: {SCHEMA} (expected deseq2|legacy)")

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

if adata.n_obs == 0 or adata.n_vars == 0:
    fail(f"empty matrix: cells={adata.n_obs} genes={adata.n_vars}")

genes = pd.Index(adata.var_names).astype(str)
if not genes.is_unique:
    duplicates = genes[genes.duplicated(keep=False)].unique().tolist()[:5]
    fail(f"gene names are not unique ({int(genes.duplicated().sum())} duplicated entries), e.g.: {duplicates}")

if sp.issparse(matrix):
    matrix = matrix.tocsr()
else:
    matrix = sp.csr_matrix(np.asarray(matrix))

# 校验覆盖全部存储条目（CSR 的 .data 即全部非零元，可能含显式零）。
data = matrix.data
if data.size:
    finite = np.isfinite(data)
    if not finite.all():
        bad = int((~finite).sum())
        fail(f"matrix has {bad} non-finite (NaN/Inf) entries; counts layer must be finite")
    negative = data < 0
    if negative.any():
        bad = int(negative.sum())
        fail(f"matrix has {bad} negative entries; counts cannot be negative (normalised layers must be re-derived, not aggregated)")
    non_integer = data != np.rint(data)
    is_integer = not non_integer.any()
    n_non_integer = int(non_integer.sum())
else:
    is_integer = True
    n_non_integer = 0

if is_integer:
    if data.size and float(np.max(np.abs(data))) > 2.0**62:
        fail("integer values exceed int64-safe range; refusing lossy aggregation")
    log_lines.append("validation\tinteger=ok\tentries_checked=%d" % data.size)
    matrix = matrix.astype(np.int64)
elif ALLOW_NON_INTEGER:
    log_lines.append(
        f"WARN matrix has {n_non_integer}/{data.size} non-integer entries; aggregating as-is "
        "(allow_non_integer=true). If values were normalised downstream must NOT model them "
        "as raw counts (plan §15.3)."
    )
    matrix = matrix.astype(np.float64)
else:
    fail(
        f"matrix has {n_non_integer} non-integer entries (of {data.size} stored); "
        "this looks like a normalised layer, not raw counts. Set allow_non_integer=true "
        "to aggregate as-is with a WARN (plan §15.3: downstream must not model it as counts)."
    )

# obs 标签缺失（NaN/None/空白串）：显式剔除进 log，缺失率 100% 才报错。
def label_mask(series):
    return series.isna().to_numpy() | (series.astype(str).str.strip() == "").to_numpy()


sample_series = adata.obs[SAMPLE_COL]
celltype_series = adata.obs[CELLTYPE_COL]
sample_missing = label_mask(sample_series)
celltype_missing = label_mask(celltype_series)
usable = ~(sample_missing | celltype_missing)
if sample_missing.all() or celltype_missing.all():
    fail(f"all cells lack a usable label (sample missing={int(sample_missing.all())}, celltype missing={int(celltype_missing.all())})")
if int((~usable).sum()):
    log_lines.append(
        f"dropped_missing\tcells={int(sample_missing.sum())}+{int(celltype_missing.sum())}\t"
        f"excluded={int((~usable).sum())}\t(unusable {SAMPLE_COL}/{CELLTYPE_COL} labels)"
    )

samples = sample_series.astype(str).to_numpy()[usable]
celltypes = celltype_series.astype(str).to_numpy()[usable]
# 标签里的 tab/换行会静默毁掉 TSV 列结构，直接拒绝。
for label_array, column_name in ((samples, SAMPLE_COL), (celltypes, CELLTYPE_COL)):
    broken = [value for value in np.unique(label_array) if ("\t" in value or "\n" in value or "\r" in value)]
    if broken:
        fail(f"obs[{column_name}] contains tab/newline characters that would corrupt the TSV output, e.g.: {broken[:3]}")

# 结构化组键：按 (sample, celltype) 元组分组；只有 ≥min_cells 的组保留。
tuple_index = {}
for sample_value, celltype_value in zip(samples, celltypes):
    tuple_index.setdefault((sample_value, celltype_value), 0)
    tuple_index[(sample_value, celltype_value)] += 1
sizes = {key: count for key, count in tuple_index.items() if count >= MIN_CELLS}
if not sizes:
    fail(f"no pseudobulk group reaches min_cells={MIN_CELLS}")
kept = sorted(sizes)
for key in sorted(set(tuple_index) - set(sizes)):
    log_lines.append(f"dropped\t{key[0]}__{key[1]}\tn_cells={tuple_index[key]}\t<min_cells={MIN_CELLS}")

# 输出 ID：默认 "<sample>__<celltype>" 拼接；不同元组拼出同一串时改百分号编码
# 回退。注意 urllib.parse.quote 对 "_"（RFC 3986 unreserved 字符）永远不转义，
# safe="" 也不行，所以必须显式 replace 把 "_"→"%5F"；这样两个分量都不含 "_"，
# "__" 只可能是连接符，拼接严格单射（fix-review R01）。
def _encode_part(value: str) -> str:
    return quote(value, safe="").replace("_", "%5F")


ids = [f"{sample_value}__{celltype_value}" for sample_value, celltype_value in kept]
id_encoding = "plain"
if len(set(ids)) != len(ids):
    ids = [_encode_part(sample_value) + "__" + _encode_part(celltype_value) for sample_value, celltype_value in kept]
    assert len(set(ids)) == len(ids), "percent-encoded ids must be injective"
    id_encoding = "percent"
    log_lines.append(
        "WARN pseudobulk id strings collide across distinct (sample, celltype) tuples; "
        "switched to percent-encoded ids (plan §15.3 audit F07; fix-review R01)"
    )
log_lines.append(f"id_encoding\t{id_encoding}")

group_of = {key: i for i, key in enumerate(kept)}
# 只有落入保留组的细胞参与聚合（缺标签、min_cells 以下的组都已排除）。
participating = np.fromiter(
    (key in group_of for key in zip(samples, celltypes)),
    dtype=bool, count=len(samples),
)
samples_used = samples[participating]
celltypes_used = celltypes[participating]
codes = np.fromiter(
    (group_of[(sample_value, celltype_value)] for sample_value, celltype_value in zip(samples_used, celltypes_used)),
    dtype=np.int64, count=len(samples_used),
)
participating_rows = np.flatnonzero(usable)[participating]

if SCHEMA == "deseq2":
    if not CONDITION_COL:
        fail("schema=deseq2 requires condition_column (the design variable taken from obs)")
    if CONDITION_COL not in adata.obs.columns:
        fail(f"obs column missing: {CONDITION_COL} (condition_column)")

    condition_series = adata.obs[CONDITION_COL]
    condition_missing = label_mask(condition_series)
    if condition_missing[participating_rows].any():
        missing_n = int(condition_missing[participating_rows].sum())
        fail(f"obs[{CONDITION_COL}] has {missing_n} missing values among aggregated cells; a design variable cannot be missing")
    condition_values = condition_series.astype(str).to_numpy()[participating_rows]
    condition_sets = [set() for _ in kept]
    for group_code, condition_value in zip(codes.tolist(), condition_values):
        condition_sets[group_code].add(condition_value)
    mixed = [kept[i] for i in range(len(kept)) if len(condition_sets[i]) > 1]
    if mixed:
        fail(
            "condition_column values are not constant within pseudobulk groups "
            f"(cells of one donor x celltype must share one condition), e.g.: {mixed[:3]}"
        )
    conditions = [sorted(condition_sets[i])[0] for i in range(len(kept))]
    log_lines.append(f"condition_column\t{CONDITION_COL}")

# 聚合保持稀疏：one-hot(cells→groups) @ 矩阵 → groups×genes，全程不稠密化。
if len(participating_rows) == adata.n_obs:
    matrix_used = matrix
else:
    matrix_used = matrix[participating_rows, :]
onehot = sp.csr_matrix(
    (np.ones(len(codes), dtype=matrix_used.dtype), (codes, np.arange(len(codes), dtype=np.int64))),
    shape=(len(kept), matrix_used.shape[0]),
)
counts_sparse = (onehot @ matrix_used).T.tocsr()  # genes × groups

# 导出按基因分块稠密化：峰值内存 ≈ 块行数 × 组数，而非 基因数 × 全量稠密矩阵。
first_col = "gene_id" if SCHEMA == "deseq2" else "gene"
with open(OUT_COUNTS, "w", encoding="utf-8", newline="\n") as handle:
    handle.write(first_col + "\t" + "\t".join(ids) + "\n")
    rows_per_chunk = max(1, min(len(genes), 10_000_000 // max(1, len(kept))))
    for start in range(0, len(genes), rows_per_chunk):
        stop = min(start + rows_per_chunk, len(genes))
        block = counts_sparse[start:stop].toarray()
        for row_offset in range(block.shape[0]):
            gene_name = genes[start + row_offset]
            values = "\t".join(str(value) for value in block[row_offset].tolist())
            handle.write(f"{gene_name}\t{values}\n")

sample_rows = []
for position, key in enumerate(kept):
    if SCHEMA == "deseq2":
        sample_rows.append((ids[position], key[0], key[1], sizes[key], conditions[position]))
    else:
        sample_rows.append((ids[position], key[0], key[1], sizes[key]))
metadata_columns = (
    ["sample_id", "sample", "celltype", "n_cells", "condition"]
    if SCHEMA == "deseq2"
    else ["pseudobulk_id", "sample", "celltype", "n_cells"]
)
pd.DataFrame(sample_rows, columns=metadata_columns).to_csv(OUT_SAMPLES, sep="\t", index=False)

log_lines.append(f"schema\t{SCHEMA}")
log_lines.append(f"groups_kept={len(kept)}\tgenes={len(genes)}\tinteger={is_integer or 'warned'}")
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write("\n".join(log_lines) + "\n")

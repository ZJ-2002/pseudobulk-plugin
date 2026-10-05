# 中文注释：审计 F07 的回归锁定套件（2026-10-05）。
# 运行：~/.cache/fix-venvs/bulk/bin/pytest tests/ -v
# 覆盖：元组组键防 "__" 碰撞（审计用例原样）、全量非零条目校验（含第 10000
# 个之后的非整数）、负值/非有限值、基因名唯一性、deseq2/legacy 表头与
# condition 列三种失败、min_cells 与缺失标签剔除、稀疏/稠密一致性与金标准对拍。
import os
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "pseudobulk_counts.py"


def run_node(tmp_path, adata, *, sample_col="donor", celltype_col="celltype",
             layer="", schema="deseq2", allow_non_integer=False,
             condition_column="", min_cells=5):
    tmp_path.mkdir(parents=True, exist_ok=True)
    input_path = tmp_path / "input.h5ad"
    adata.write_h5ad(input_path)
    counts = tmp_path / "counts.tsv"
    samples = tmp_path / "samples.tsv"
    log = tmp_path / "pseudobulk.log"
    env = dict(os.environ)
    env.update({
        "PB_SAMPLE_COL": sample_col,
        "PB_CELLTYPE_COL": celltype_col,
        "PB_MIN_CELLS": str(min_cells),
        "PB_LAYER": layer,
        "PB_SCHEMA": schema,
        "PB_ALLOW_NON_INTEGER": "true" if allow_non_integer else "false",
        "PB_CONDITION_COLUMN": condition_column,
        "AUTONOMICS_INPUT0": str(input_path),
        "AUTONOMICS_OUTPUT0": str(counts),
        "AUTONOMICS_OUTPUT1": str(samples),
        "AUTONOMICS_OUTPUT2": str(log),
    })
    proc = subprocess.run([sys.executable, str(SCRIPT)], env=env,
                          capture_output=True, text=True)
    return proc, counts, samples, log


def make_adata(sample_labels, celltype_labels, X, *, gene_names=None,
               condition_labels=None):
    # 注：anndata 0.13 + pandas 3 写 object/categorical obs 列会丢 codes，
    # fixture 一律用 pandas 默认 StringDtype 纯串（与引擎 Rust 写的纯字符串
    # dataset 回读路径一致）；缺失标签用空串或 float NaN，不用 None。
    obs = pd.DataFrame({
        "donor": list(sample_labels),
        "celltype": list(celltype_labels),
    }, index=[f"cell{i}" for i in range(len(sample_labels))])
    if condition_labels is not None:
        obs["condition"] = list(condition_labels)
    var = pd.DataFrame(index=list(gene_names) if gene_names is not None
                       else [f"g{j}" for j in range(X.shape[1])])
    return ad.AnnData(X=X, obs=obs, var=var)


def sparse_counts(n_cells, n_genes, seed=7):
    rng = np.random.default_rng(seed)
    return sp.csr_matrix(rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.int64))


def read_counts(path):
    table = pd.read_csv(path, sep="\t", index_col=0)
    return table


def read_samples(path):
    return pd.read_csv(path, sep="\t")


def test_tuple_keys_survive_separator_collision(tmp_path):
    """审计原始用例：("donor__state","type") 与 ("donor","state__type") 不得合并。"""
    n = 12
    donors = ["donor__state"] * n + ["donor"] * n
    celltypes = ["type"] * n + ["state__type"] * n
    X = sp.csr_matrix((2 * n, 2), dtype=np.int64)
    X[:n, 0] = 3  # 第一组在 g0 上累计
    X[n:, 1] = 4  # 第二组在 g1 上累计
    adata = make_adata(donors, celltypes, X,
                       condition_labels=["ctrl"] * n + ["stim"] * n)
    proc, counts_path, samples_path, log_path = run_node(
        tmp_path, adata, condition_column="condition")
    assert proc.returncode == 0, proc.stderr

    counts = read_counts(counts_path)
    samples = read_samples(samples_path)
    assert counts.shape == (2, 2)
    # 两组各 12 细胞：若被拼接键合并会只剩一组 24 细胞。
    assert sorted(samples["n_cells"].tolist()) == [12, 12]
    assert set(map(tuple, samples[["sample", "celltype"]].values)) == {
        ("donor__state", "type"), ("donor", "state__type")}
    # fix-review R01：必须检查原始 TSV 表头——pandas 会把重复列名自动重命名
    # 成 ".1"，经 pandas 读回的 counts.columns 看不到真实重复。
    raw_header = counts_path.read_text().splitlines()[0].split("\t")
    raw_ids = [name for name in raw_header[1:]]
    assert len(raw_ids) == 2 and len(set(raw_ids)) == 2, f"raw header ids not unique: {raw_ids}"
    assert all(".1" not in name for name in counts.columns), "pandas renamed duplicate columns"
    # counts 原始表头与 metadata 的 sample_id 逐字一致（DESeq2 身份契约）。
    assert samples["sample_id"].tolist() == raw_ids
    assert counts.columns.tolist() == raw_ids
    # "_" 编码为 %5F 后，连接符 "__" 不再有歧义。
    assert sorted(raw_ids) == sorted(["donor%5F%5Fstate__type", "donor__state%5F%5Ftype"])
    # kept 排序：("donor","state__type") 在 ("donor__state","type") 之前。
    assert counts.loc["g0"].tolist() == [0, 3 * n]
    assert counts.loc["g1"].tolist() == [4 * n, 0]
    assert "id_encoding\tpercent" in log_path.read_text()


def test_plain_ids_keep_readable_join(tmp_path):
    donors = ["d1", "d1", "d2", "d2"] * 5
    celltypes = ["A", "B", "A", "B"] * 5
    adata = make_adata(donors, celltypes, sparse_counts(20, 6),
                       condition_labels=["ctrl", "ctrl", "stim", "stim"] * 5)
    proc, counts_path, samples_path, log_path = run_node(
        tmp_path, adata, condition_column="condition")
    assert proc.returncode == 0, proc.stderr
    counts = read_counts(counts_path)
    assert counts.columns.tolist() == ["d1__A", "d1__B", "d2__A", "d2__B"]
    assert "id_encoding\tplain" in log_path.read_text()


def test_non_integer_after_first_10000_entries_hard_error(tmp_path):
    """稀疏校验必须覆盖第 10000 个存储条目之后的非整数。"""
    n_cells, n_genes = 2000, 40  # 2000×40 随机 Poisson ≫ 10000 个非零
    rng = np.random.default_rng(11)
    dense = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float64)
    donors = [f"d{i % 4}" for i in range(n_cells)]
    celltypes = [f"t{i % 2}" for i in range(n_cells)]
    conditions = ["ctrl", "stim"] * (n_cells // 2)
    X = sp.csr_matrix(dense)
    assert X.nnz > 10000
    X.data[15000] = 0.5  # 埋在老代码 10000 探针之外
    adata = make_adata(donors, celltypes, X, condition_labels=conditions)

    proc, *_ = run_node(tmp_path, adata, condition_column="condition")
    assert proc.returncode == 2
    assert "non-integer entries" in proc.stderr

    proc2, counts2, _, log2 = run_node(tmp_path / "allow", adata,
                                       condition_column="condition",
                                       allow_non_integer=True)
    assert proc2.returncode == 0, proc2.stderr
    assert "non-integer entries" in log2.read_text()
    assert "allow_non_integer=true" in log2.read_text()


def test_negative_and_nonfinite_are_hard_errors(tmp_path):
    donors = ["d1"] * 8 + ["d2"] * 8
    celltypes = ["A"] * 16
    conditions = ["ctrl"] * 16

    X = sparse_counts(16, 5)
    X.data[3] = -2.0
    proc, *_ = run_node(tmp_path, make_adata(donors, celltypes, X, condition_labels=conditions),
                        condition_column="condition")
    assert proc.returncode == 2
    assert "negative" in proc.stderr

    X2 = sparse_counts(16, 5, seed=3).astype(np.float64)
    X2.data[2] = np.nan
    proc2, *_ = run_node(tmp_path / "nan",
                         make_adata(donors, celltypes, X2, condition_labels=conditions),
                         condition_column="condition")
    assert proc2.returncode == 2
    assert "non-finite" in proc2.stderr


def test_duplicate_gene_names_error(tmp_path):
    adata = make_adata(["d1"] * 10, ["A"] * 10, sparse_counts(10, 4),
                       gene_names=["g1", "g2", "g2", "g3"],
                       condition_labels=["ctrl"] * 10)
    proc, *_ = run_node(tmp_path, adata, condition_column="condition")
    assert proc.returncode == 2
    assert "not unique" in proc.stderr


def test_deseq2_schema_headers_and_condition(tmp_path):
    donors = ["d1", "d1", "d2", "d2"] * 5
    celltypes = ["A", "B", "A", "B"] * 5
    conditions = ["ctrl", "ctrl", "stim", "stim"] * 5
    adata = make_adata(donors, celltypes, sparse_counts(20, 6),
                       condition_labels=conditions)
    proc, counts_path, samples_path, _ = run_node(tmp_path, adata,
                                                  condition_column="condition")
    assert proc.returncode == 0, proc.stderr
    header = counts_path.read_text().splitlines()[0].split("\t")
    assert header[0] == "gene_id"
    samples = read_samples(samples_path)
    counts_df = read_counts(counts_path)
    assert samples.columns.tolist() == ["sample_id", "sample", "celltype", "n_cells", "condition"]
    assert samples["sample_id"].tolist() == counts_df.columns.tolist()
    assert dict(zip(samples["sample_id"], samples["condition"])) == {
        "d1__A": "ctrl", "d1__B": "ctrl", "d2__A": "stim", "d2__B": "stim"}


def test_deseq2_requires_condition_column_param(tmp_path):
    adata = make_adata(["d1"] * 10, ["A"] * 10, sparse_counts(10, 4))
    proc, *_ = run_node(tmp_path, adata, schema="deseq2")
    assert proc.returncode == 2
    assert "condition_column" in proc.stderr


def test_deseq2_condition_column_must_exist_in_obs(tmp_path):
    adata = make_adata(["d1"] * 10, ["A"] * 10, sparse_counts(10, 4))
    proc, *_ = run_node(tmp_path, adata, condition_column="treatment")
    assert proc.returncode == 2
    assert "obs column missing: treatment" in proc.stderr


def test_deseq2_condition_must_be_constant_within_group(tmp_path):
    donors = ["d1"] * 10
    celltypes = ["A"] * 10
    conditions = ["ctrl"] * 5 + ["stim"] * 5  # 同组内两种 condition = 数据矛盾
    adata = make_adata(donors, celltypes, sparse_counts(10, 4),
                       condition_labels=conditions)
    proc, *_ = run_node(tmp_path, adata, condition_column="condition")
    assert proc.returncode == 2
    assert "not constant" in proc.stderr


def test_legacy_schema_keeps_old_headers(tmp_path):
    donors = ["d1", "d2"] * 5
    celltypes = ["A"] * 10
    adata = make_adata(donors, celltypes, sparse_counts(10, 4))  # 不给 condition
    proc, counts_path, samples_path, _ = run_node(tmp_path, adata, schema="legacy")
    assert proc.returncode == 0, proc.stderr
    header = counts_path.read_text().splitlines()[0].split("\t")
    assert header[0] == "gene"
    samples = read_samples(samples_path)
    assert samples.columns.tolist() == ["pseudobulk_id", "sample", "celltype", "n_cells"]
    assert samples["pseudobulk_id"].tolist() == ["d1__A", "d2__A"]


def test_min_cells_drop_and_missing_labels_excluded(tmp_path):
    donors = (["d1"] * 8 + ["d2"] * 8 + ["d3"] * 3          # d3 整组 < min_cells
              + ["d1"] * 6)                                  # 缺 celltype 的 6 个
    celltypes = ["A"] * 8 + ["A"] * 8 + ["A"] * 3 + [""] * 6  # 空串 = 标签缺失
    conditions = ["ctrl"] * 25  # 缺标签细胞已被剔除，其 condition 值不参与校验
    X = sparse_counts(25, 5)
    X[:8, 0] = 5  # d1__A 的金标准：8 细胞 × 5
    adata = make_adata(donors, celltypes, X, condition_labels=conditions)
    proc, counts_path, samples_path, log_path = run_node(tmp_path, adata,
                                                         condition_column="condition")
    assert proc.returncode == 0, proc.stderr
    log_text = log_path.read_text()
    assert "dropped\td3__A\tn_cells=3\t<min_cells=5" in log_text
    assert "dropped_missing" in log_text
    samples = read_samples(samples_path)
    assert samples["sample"].tolist() == ["d1", "d2"]
    assert samples["n_cells"].tolist() == [8, 8]
    counts = read_counts(counts_path)
    assert counts.loc["g0", "d1__A"] == 8 * 5  # 缺标签细胞不并入任何组


def test_nan_donor_labels_excluded(tmp_path):
    donors = [1.0] * 8 + [np.nan] * 4  # float 列的 NaN 走 isna 缺失路径
    celltypes = ["A"] * 12
    adata = make_adata(donors, celltypes, sparse_counts(12, 4))
    proc, counts_path, samples_path, log_path = run_node(tmp_path, adata, schema="legacy")
    assert proc.returncode == 0, proc.stderr
    samples = read_samples(samples_path)
    # 样本表里 sample 列是字符串 "1.0"；read_csv 数值推断后需转回比较。
    assert samples["sample"].astype(str).tolist() == ["1.0"]
    assert samples["n_cells"].tolist() == [8]
    assert "dropped_missing" in log_path.read_text()


def test_sparse_dense_agree_and_match_manual_groupby(tmp_path):
    rng = np.random.default_rng(23)
    n = 60
    donors = [f"d{i % 3}" for i in range(n)]
    celltypes = [f"t{i % 2}" for i in range(n)]
    dense = rng.poisson(3.0, size=(n, 8)).astype(np.int64)

    adata_sparse = make_adata(donors, celltypes, sp.csr_matrix(dense))
    adata_dense = make_adata(donors, celltypes, dense.copy())
    proc_a, counts_a, samples_a, _ = run_node(tmp_path, adata_sparse, schema="legacy")
    proc_b, counts_b, _, _ = run_node(tmp_path / "dense", adata_dense, schema="legacy")
    assert proc_a.returncode == 0, proc_a.stderr
    assert proc_b.returncode == 0, proc_b.stderr

    counts_a_df = read_counts(counts_a)
    counts_b_df = read_counts(counts_b)
    pd.testing.assert_frame_equal(counts_a_df, counts_b_df)

    # 金标准：独立用 pandas groupby 重算（与脚本 one-hot 路径无共享代码）。
    frame = pd.DataFrame({"donor": donors, "celltype": celltypes})
    for gene_index, gene_name in enumerate(counts_a_df.index):
        frame["value"] = dense[:, gene_index]
        golden = frame.groupby(["donor", "celltype"])["value"].sum()
        for column in counts_a_df.columns:
            sample_value, celltype_value = column.split("__", 1)
            assert counts_a_df.loc[gene_name, column] == golden.loc[(sample_value, celltype_value)]


def test_named_layer_used_when_requested(tmp_path):
    donors = ["d1"] * 6
    celltypes = ["A"] * 6
    conditions = ["ctrl"] * 6
    X = sparse_counts(6, 4, seed=1)
    counts_layer = sparse_counts(6, 4, seed=2) * 10
    adata = make_adata(donors, celltypes, X, condition_labels=conditions)
    adata.layers["counts"] = counts_layer
    proc_raw, counts_raw, _, _ = run_node(tmp_path, adata, condition_column="condition")
    proc_layer, counts_layer_path, _, log_path = run_node(tmp_path / "layer", adata,
                                                          condition_column="condition",
                                                          layer="counts")
    assert proc_raw.returncode == 0 and proc_layer.returncode == 0
    assert "layer\tcounts" in log_path.read_text()
    layer_df = read_counts(counts_layer_path)
    assert layer_df.values.sum() == int(counts_layer.sum())
    frame = pd.DataFrame({
        "donor": donors, "celltype": celltypes,
        "v": counts_layer[:, 0].toarray().ravel(),  # 只对拍 g0 这一列
    })
    golden = frame.groupby(["donor", "celltype"])["v"].sum()
    assert layer_df.loc["g0", "d1__A"] == golden.loc[("d1", "A")]

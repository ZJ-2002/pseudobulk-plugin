# 中文注释：F07 资源基准——20k 细胞 × 500 基因量级的 pseudobulk 聚合。
# 运行：~/.cache/fix-venvs/bulk/bin/python tests/bench_pseudobulk.py
# 输出：墙钟时间与子进程峰值 RSS（resource.ru_maxrss，KB，取本进程所有
# 子进程的最大值；基准脚本只跑一个节点子进程，即节点峰值）。
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "pseudobulk_counts.py"

N_CELLS = 20_000
N_GENES = 500
N_DONORS = 5
N_CELLTYPES = 4  # 20 组 × 1000 细胞

rng = np.random.default_rng(42)
X = sp.csr_matrix(rng.poisson(1.5, size=(N_CELLS, N_GENES)).astype(np.int64))
donors = np.repeat([f"d{i}" for i in range(N_DONORS)], N_CELLS // N_DONORS)
celltypes = np.tile(np.repeat([f"t{j}" for j in range(N_CELLTYPES)], N_CELLS // (N_DONORS * N_CELLTYPES)), N_DONORS)
obs = pd.DataFrame({
    "donor": list(map(str, donors)),
    "celltype": list(map(str, celltypes)),
    # condition 挂在供者上（组内恒定），偶数供者对照、奇数供者处理。
    "condition": ["ctrl" if int(d[1:]) % 2 == 0 else "stim" for d in donors],
}, index=[f"cell{i}" for i in range(N_CELLS)])
adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"gene{i}" for i in range(N_GENES)]))

with tempfile.TemporaryDirectory(prefix="pb_bench_") as work:
    work = Path(work)
    input_path = work / "input.h5ad"
    adata.write_h5ad(input_path)
    env = dict(os.environ)
    env.update({
        "PB_SAMPLE_COL": "donor", "PB_CELLTYPE_COL": "celltype",
        "PB_MIN_CELLS": "10", "PB_LAYER": "",
        "PB_SCHEMA": "deseq2", "PB_ALLOW_NON_INTEGER": "false",
        "PB_CONDITION_COLUMN": "condition",
        "AUTONOMICS_INPUT0": str(input_path),
        "AUTONOMICS_OUTPUT0": str(work / "counts.tsv"),
        "AUTONOMICS_OUTPUT1": str(work / "samples.tsv"),
        "AUTONOMICS_OUTPUT2": str(work / "pseudobulk.log"),
    })
    start = time.perf_counter()
    proc = subprocess.run([sys.executable, str(SCRIPT)], env=env, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    peak_rss_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if proc.returncode != 0:
        print(proc.stderr)
        sys.exit(1)
    n_groups = len((work / "samples.tsv").read_text().splitlines()) - 1

print(f"input: {N_CELLS} cells x {N_GENES} genes, nnz={X.nnz}, groups={n_groups}")
print(f"wall time: {elapsed:.2f} s")
print(f"child peak RSS: {peak_rss_kb / 1024:.1f} MiB")

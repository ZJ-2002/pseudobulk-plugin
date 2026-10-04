# pseudobulk

单细胞 h5ad → 供者×细胞类型假体块计数。补齐主计划 §15.3 的缺口：
差异分析的统计单位是**供者**不是细胞；引擎原 registry 无 pseudobulk
聚合工厂。2026-10-04 增补 `h5ad_export_matrix`：X 层长表导出（引擎
原生 h5ad 工厂只有 obs 层，X 无导出路径的补位）。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `pseudobulk_counts` | h5ad（AnnData） | `pseudobulk_counts.tsv`（基因 × `<sample>__<celltype>`）、`pseudobulk_samples.tsv`（含 `n_cells`）、log |
| `h5ad_export_matrix` | h5ad（任意来源） | `x_matrix.long.parquet`（cell_id/gene/value）、`obs.tsv`、`var.tsv`、log |

参数：`pseudobulk_counts` 的 `sample_col`、`celltype_col`（必填），
`min_cells`（默认 10），`layer`（可选，默认 `.X`），`id_schema`
（默认 `{sample}__{celltype}`）；`h5ad_export_matrix` 的 `layer`
（默认 `X`）、`chunk`（默认 2000 细胞/行组）。

## 布局

```text
pseudobulk/
├── manifest.toml
├── scripts/pseudobulk_counts.py
├── scripts/h5ad_export_matrix.py
├── Dockerfile
└── README.md
```

镜像：`localhost/autonomics/biotools-py@sha256:21e14c1582a9d4c261c0a23a11864293313ee2b49b9b7d85bcabc89210b48c7e`。

## 口径（对齐 §15.3/§16）

- one-hot 稀疏矩阵乘聚合（`onehot @ matrix`），不稠密化。
- 计数非整数 → **WARN 不报错**：归一化矩阵被聚合时结果仍是合法
  均值矩阵，但下游 DESeq2 只能喂整数计数层；边界写进 log。
- `n_cells < min_cells` 的组剔除并逐组记 log（低细胞数组不进统计，
  防单细胞漂移伪装成生物学）。
- 输出样本表带 `n_cells` 列，供下游 LMM/DESeq2 检查每供者贡献。

## 冒烟记录

- `pseudobulk_counts`：6 供者 Poisson 梯度合成数据 D1≈214 → D6≈456
  单调恢复；稀疏/稠密矩阵路径一致。
- `h5ad_export_matrix`（2026-10-04，biotools-py 容器内实跑）：
  引擎自产 lzf-CSR h5ad（sham_reannot，10794×21355）→ 24,376,851
  非零元、cell_id/基因符号正确；宿主 anndata 0.13 写的 h5ad（10×5）
  → 50 条全量导出。双向可读（lzf filter 32000 = h5py 官方插件 ID）。

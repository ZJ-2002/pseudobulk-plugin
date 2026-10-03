# pseudobulk

单细胞 h5ad → 供者×细胞类型假体块计数。补齐主计划 §15.3 的缺口：
差异分析的统计单位是**供者**不是细胞；引擎原 registry 无 pseudobulk
聚合工厂。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `pseudobulk_counts` | h5ad（AnnData） | `pseudobulk_counts.tsv`（基因 × `<sample>__<celltype>`）、`pseudobulk_samples.tsv`（含 `n_cells`）、log |

参数：`sample_col`、`celltype_col`（必填），`min_cells`（默认 10），
`layer`（可选，默认 `.X`），`id_schema`（默认 `{sample}__{celltype}`）。

## 布局

```text
pseudobulk/
├── manifest.toml
├── scripts/pseudobulk_counts.py
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

6 供者 Poisson 梯度合成数据：D1≈214 → D6≈456 单调恢复；稀疏/稠密
矩阵路径一致。

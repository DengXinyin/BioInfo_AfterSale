# 差异丰度可视化

对应 `R/04_differential.R`。

- `metage_diff_boxplot()`：单特征箱线图；已有 P 值通过 `p_value` 传入，否则可明确选择 Wilcoxon/ANOVA。
- `metage_stamp_plot()`：读取已有 STAMP 差异、置信区间和 P 值列，不隐式改换 raw/adjusted P。
- `metage_random_forest_importance()`：读取已有变量重要性表，不重新拟合模型。
- `metage_lda_score_plot()`：读取 LEfSe marker 表，`cutoff` 可配置。
- `metage_significant_heatmap()`：按已有结果筛选显著特征后画热图。

```r
p <- metage_lda_score_plot(marker, feature = "Taxon", group = "Class",
                            score = "LDA", cutoff = 3)
```

镜像内 R 脚本只包含 LDA 条形图，不包含 cladogram；需要 cladogram 时使用仓库根包的 `LEFSE_cladogram_plot()`。

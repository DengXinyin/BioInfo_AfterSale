# Alpha 多样性与排序分析

对应 `R/03_alpha_ordination.R`。

`metage_alpha_indices()` 从原始计数矩阵计算 Observed、Chao1、ACE、Shannon 与 Gini-Simpson。`metage_alpha_boxplot()` 可执行双组 Wilcoxon 或多组 ANOVA，并同时返回图、P 值和绘图数据。

`metage_ordination()` 支持 PCA、PCoA、NMDS。返修已有结果时建议传 `coordinates`，这样只重绘、不重算；只有没有既有坐标时才传 `abundance`。分组向量建议以样本名命名。小样本椭圆优先使用 `ggforce::geom_mark_ellipse()`。

```r
ord <- metage_ordination(abundance, group, method = "PCoA", ellipse = TRUE)
metage_save_plot(ord$plot, "Output/pcoa.pdf", 8, 6)
```

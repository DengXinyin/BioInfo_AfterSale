# 点阵相关热图教程

本示例使用虚构的“样本 × 特征”数据，展示特征与环境变量的 Spearman 相关性。行是样本，列是变量，行名是样本标识；示例只输出 PDF。

```r
library(BioInfoAfterSale)
set.seed(20260903)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
feature <- as.data.frame(matrix(rnorm(12 * 16, 8, 1.5), 12, 16, dimnames = list(sample_id, paste0("Feature_", 1:16))))
property <- as.data.frame(matrix(rnorm(12 * 8), 12, 8, dimnames = list(sample_id, paste0("Property_", 1:8))))
style <- choose_plot_style(font_family = "Times New Roman", theme = "classic", dpi = 300, figure_width = 11, figure_height = 8, title = list(size = 18, bold = TRUE))
dir.create("docs/images", showWarnings = FALSE)
p <- plot_dot_heatmap(feature, property, method = "spearman", title = "Feature vs property", style = style, output_file = "docs/images/Heatmap_dot.pdf")
p
```

结果：[Heatmap_dot.pdf](images/Heatmap_dot.pdf)

点的颜色表示相关方向和强度，点的面积表示 `|r|`。至少需要 3 个共同样本；常数列或有效值不足的组合显示为空值。

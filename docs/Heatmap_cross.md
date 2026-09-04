# 矩形交叉相关热图教程

本示例使用虚构的物种丰度表和环境变量表。两张表均为“样本在行、变量在列”，通过相同的样本行名配对；示例只输出 PDF。

```r
library(BioInfoAfterSale)
set.seed(20260903)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
species <- as.data.frame(matrix(abs(rnorm(12 * 24, 20, 6)), 12, 24, dimnames = list(sample_id, paste0("Species_", 1:24))))
environment <- as.data.frame(matrix(rnorm(12 * 9), 12, 9, dimnames = list(sample_id, paste0("Property_", 1:9))))
style <- choose_plot_style(font_family = "Times New Roman", theme = "classic", dpi = 300, figure_width = 11, figure_height = 9, title = list(size = 18, bold = TRUE))
dir.create("docs/images", showWarnings = FALSE)
p <- plot_cross_heatmap(species, environment, method = "spearman", title = "Species vs properties", style = style, output_file = "docs/images/Heatmap_cross.pdf")
p
```

结果：[Heatmap_cross.pdf](images/Heatmap_cross.pdf)

色块表示 Spearman 相关系数：蓝色为负相关，红色为正相关，灰色为无法计算的组合。

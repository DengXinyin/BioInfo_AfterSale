# 圆形相关热图教程

本示例使用虚构的物种丰度和环境变量数据，并采用 `circlize::circos.heatmap()` 绘制同心环热图。每个扇区对应一个物种，每个同心环对应一个环境变量；数据通过样本行名配对，示例只输出 PDF。

```r
library(BioInfoAfterSale)
set.seed(20260903)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
species <- as.data.frame(matrix(abs(rnorm(12 * 18, 15, 5)), 12, 18, dimnames = list(sample_id, paste0("Species_", 1:18))))
environment <- as.data.frame(matrix(rnorm(12 * 7), 12, 7, dimnames = list(sample_id, paste0("Property_", 1:7))))
style <- choose_plot_style(font_family = "Times New Roman", theme = "classic", dpi = 300, figure_width = 12, figure_height = 12, title = list(size = 18, bold = TRUE), legend_title = list(size = 14), legend_text = list(size = 12))
dir.create("docs/images", showWarnings = FALSE)
rho <- plot_circular_heatmap(species, environment, method = "spearman", title = "Circular species-property correlation", style = style, output_file = "docs/images/Heatmap_circular.pdf")
rho[1:3, 1:3]
```

结果：[Heatmap_circular.pdf](images/Heatmap_circular.pdf)

圆形图保留平均连接、欧氏距离聚类和内部树状图；颜色使用 `-1` 到 `1` 的连续相关系数标尺。绘图结束后会自动清理圆形绘图状态。

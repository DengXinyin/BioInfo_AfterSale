# 矩形交叉相关热图

使用色块展示两个变量集合之间的 Spearman 相关系数。

```r
library(BioInfoAfterSale)
set.seed(20260902)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
make_table <- function(prefix, n) as.data.frame(matrix(rnorm(12 * n), 12, n, dimnames = list(sample_id, paste0(prefix, 1:n))))
style <- choose_plot_style(font_family = "sans", theme = "classic", dpi = 300, figure_width = 9, figure_height = 7)
p <- plot_cross_heatmap(make_table("Taxon_", 10), make_table("Property_", 5), title = "Synthetic taxa vs soil properties", style = style, output_file = "docs/images/plot-cross-heatmap.png")
p
```

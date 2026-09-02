# 圆形交叉相关热图

以极坐标紧凑展示变量集合之间的相关矩阵，适合报告版式。

```r
library(BioInfoAfterSale)
set.seed(20260902)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
make_table <- function(prefix, n) as.data.frame(matrix(rnorm(12 * n), 12, n, dimnames = list(sample_id, paste0(prefix, 1:n))))
style <- choose_plot_style(font_family = "sans", theme = "classic", dpi = 300, figure_width = 9, figure_height = 7)
p <- plot_circular_heatmap(make_table("Taxon_", 10), make_table("Property_", 5), title = "Circular cross-correlation heatmap", style = style, output_file = "docs/images/plot-circular-heatmap.png")
p
```

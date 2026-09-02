# 点阵相关热图

适合同时展示相关方向和强度。输入数据须为“样本在行、变量在列”，并带有相同的样本行名。

```r
library(BioInfoAfterSale)
set.seed(20260902)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
make_table <- function(prefix, n) as.data.frame(matrix(rnorm(12 * n), 12, n, dimnames = list(sample_id, paste0(prefix, 1:n))))
style <- choose_plot_style(font_family = "sans", theme = "classic", dpi = 300, figure_width = 9, figure_height = 7)
p <- plot_dot_heatmap(make_table("Metabolite_", 8), make_table("Property_", 5), title = "Synthetic metabolites vs soil properties", style = style, output_file = "docs/images/plot-dot-heatmap.png")
p
```

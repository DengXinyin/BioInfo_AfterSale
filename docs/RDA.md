# RDA 排序图

`plot_rda()` 对响应矩阵执行 Hellinger 变换、VIF 筛选和置换检验，并绘制样本分组与环境变量箭头。

```r
library(BioInfoAfterSale)
set.seed(20260902)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
make_table <- function(prefix, n) as.data.frame(matrix(rnorm(12 * n), 12, n, dimnames = list(sample_id, paste0(prefix, 1:n))))
response <- abs(make_table("Taxon_", 10)); environment <- make_table("Property_", 5)
groups <- setNames(rep(c("A", "B"), each = 6), sample_id)
style <- choose_plot_style(font_family = "Times New Roman", theme = "classic", dpi = 300, figure_width = 9, figure_height = 7)
rda_result <- plot_rda(response, environment, group = groups, title = "Synthetic community constrained by soil properties", style = style, output_file = "docs/images/plot-rda.pdf")
rda_result$plot
rda_result$selected_variables
rda_result$anova
```

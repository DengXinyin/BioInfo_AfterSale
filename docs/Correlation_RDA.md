# 相关热图与 RDA 总览

本模块已拆分为 4 个独立函数脚本，并配有同名教程：

- [plot_dot_heatmap](plot_dot_heatmap.md)
- [plot_cross_heatmap](plot_cross_heatmap.md)
- [plot_circular_heatmap](plot_circular_heatmap.md)
- [plot_rda](plot_rda.md)

下面保留完整的联合示例，便于一次性复现全部结果图。

本页使用完全模拟的数据演示 3 种相关热图和 RDA。示例中的变量名、样本名和数值均为新生成内容，不来自任何项目数据。

## 统一样式与模拟数据

```r
library(BioInfoAfterSale)
set.seed(20260902)
sample_id <- paste0("Sample_", sprintf("%02d", 1:12))
make_table <- function(prefix, n) {
  z <- matrix(rnorm(12 * n), 12, n, dimnames = list(sample_id, paste0(prefix, 1:n)))
  as.data.frame(z, check.names = FALSE)
}
exposure <- make_table("Metabolite_", 8)
soil <- make_table("Property_", 5)
microbe <- abs(make_table("Taxon_", 10))
groups <- setNames(rep(c("A", "B"), each = 6), sample_id)

style <- choose_plot_style(
  font_family = "sans", theme = "classic", dpi = 300,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE),
  axis_text = list(size = 11), legend_text = list(size = 10),
  group_palette = c("#18BFC2", "#F8766D")
)
dir.create("correlation_rda_results", showWarnings = FALSE)
```

## 1. 点阵相关热图

适合变量数量较多、希望同时表达相关系数方向和绝对值大小的场景。

```r
p_dot <- plot_dot_heatmap(
  exposure, soil, title = "Synthetic metabolites vs soil properties",
  style = style, output_file = "correlation_rda_results/dot_heatmap.pdf"
)
p_dot
```

颜色表示相关系数，点的面积表示 `|r|`。函数默认按样本名取两个表的交集，并要求至少 3 个共同样本。

## 2. 矩形交叉热图

适合完整查看两个变量集合之间的相关矩阵，例如微生物分类单元与环境因子。

```r
p_cross <- plot_cross_heatmap(
  microbe, soil, title = "Synthetic taxa vs soil properties",
  style = style, output_file = "correlation_rda_results/cross_heatmap.pdf"
)
p_cross
```

## 3. 圆形交叉热图

适合变量较多、需要紧凑展示交叉矩阵的报告版式。它与矩形热图使用同一相关计算和颜色标尺。

```r
p_circular <- plot_circular_heatmap(
  microbe, soil, title = "Circular cross-correlation heatmap",
  style = style, output_file = "correlation_rda_results/circular_heatmap.pdf"
)
p_circular
```

## RDA 排序图

`plot_rda()` 要求响应表和环境表均为“样本在行、变量在列”，并且具有相同的样本行名。默认对响应数据做 Hellinger 变换，逐步移除最高 VIF 变量，并在最多保留 3 个解释变量后进行 999 次置换检验。

```r
rda_result <- plot_rda(
  response = microbe, environment = soil, group = groups,
  title = "Synthetic community constrained by soil properties",
  style = style, output_file = "correlation_rda_results/rda.pdf"
)
rda_result$plot
rda_result$selected_variables
rda_result$anova
```

若响应矩阵已经完成适当变换，可使用 `transform = "none"`。`rda_result$model`、`site_scores` 和 `environmental_vectors` 可用于后续导出或自定义标注。

# GSEA 可视化示例

使用 `plot_gsea()` 绘制运行富集分数、基因集命中位置和排序统计量三个面板。
`hits` 可以是命中位置，也可以是与排序向量等长的逻辑向量。

所有字体优先使用 `Times New Roman`。NES、P value 和 Adjusted P 只标注在最上方的运行富集分数面板。

## 示例数据与绘图

```r
library(BioInfoAfterSale)
after_sale_style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 400,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 14, bold = TRUE),
  axis_text = list(size = 12, bold = TRUE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)

set.seed(2026)
n_rank <- 12000
ranked_metric <- sort(rnorm(n_rank), decreasing = TRUE)
hit_positions <- sort(sample(seq_len(n_rank), 260))
is_hit <- seq_len(n_rank) %in% hit_positions
hit_weight <- abs(ranked_metric) * is_hit
running_score <- cumsum(hit_weight / sum(hit_weight) - (!is_hit) / sum(!is_hit))

plot_gsea(
  running_score = running_score, hits = hit_positions, ranked_metric = ranked_metric,
  title = "Peptide cross-linking",
  statistics = c(NES = 2.22, `P value` = 6.1e-09, `Adjusted P` = 2.8e-06),
  colors = c("#E25659", "#335372"), style = after_sale_style,
  output_file = "docs/images/GSEA.png", width = 9.2, height = 6.8
)
```

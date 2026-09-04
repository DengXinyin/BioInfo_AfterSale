# 火山图可视化示例

使用 `plot_volcano()` 展示差异基因的 log2 Fold Change 与显著性。

所有字体优先使用 `Times New Roman`。组别图例放在主图内部，并显示每组数量，例如 `Up(206)`。

## 示例数据与绘图

```r
library(BioInfoAfterSale)
after_sale_style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 400,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 14, bold = TRUE),
  axis_text = list(size = 12, bold = TRUE), legend_text = list(size = 12, bold = TRUE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)

set.seed(2026)
n_gene <- 5000
deg_data <- data.frame(
  GeneID = paste0("Gene", seq_len(n_gene)),
  log2FC = rnorm(n_gene, 0, 1.8),
  pvalue = pmin(runif(n_gene), 10^(-rexp(n_gene, rate = 0.8)))
)

plot_volcano(
  deg_data, log2fc = "log2FC", pvalue = "pvalue", fc_cutoff = 1,
  p_cutoff = 0.05, status_colors = c("#E25659", "#D4D4D4", "#335372"),
  title = "TA4% vs Control", legend_title = NULL, legend_inside = c(0.82, 0.82),
  style = after_sale_style,
  output_file = "docs/images/Volcano.pdf"
)
```

已有差异状态列时，将 `status = "level"` 传给函数，并保证状态值为 `Up`、`Down` 或 `Not sig`。

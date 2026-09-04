# 三集合 Venn 图可视化示例

`plot_venn()` 接收一个包含 3 个命名向量的列表，并绘制集合之间的基因重叠关系。

所有字体优先使用 `Times New Roman`。Venn 图不绘制边框和网格线，集合名直接标注在对应圆圈边上，不使用独立图例。

## 示例数据与绘图

```r
library(BioInfoAfterSale)
after_sale_style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 400,
  figure_width = 9, figure_height = 7,
  title = list(size = 20, bold = TRUE), axis_title = list(size = 18, bold = TRUE),
  axis_text = list(show = FALSE), data_label = list(size = 16, bold = TRUE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)

set.seed(2026)
gene_universe <- paste0("Gene", 1:5000)
venn_sets <- list(
  `TA4% vs TA0%` = sample(gene_universe, 1800),
  `TA4% vs Control` = sample(gene_universe, 2200),
  `TA0% vs Control` = sample(gene_universe, 1500)
)

plot_venn(venn_sets, title = "Differentially Expressed Gene Overlap",
          style = after_sale_style, output_file = "docs/images/Venn.pdf",
          width = 8, height = 7)
```

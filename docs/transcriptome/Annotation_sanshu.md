# 转录组功能注释汇总图

本教程演示注释花瓣图、功能分类柱状图和物种组成饼图。数据均为虚构汇总值；`Database` 表示数据库，`Category` 表示功能类别或物种，`Count` 表示基因数量。数值型类别按总数量降序排列，图例不带黑色外框。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 13),
  axis_text = list(size = 10), legend = list(position = "right", frame = FALSE)
)

annotation_summary <- data.frame(
  Database = c("NR", "Swiss-Prot", "GO", "KEGG", "KOG", "Pfam", "TF"),
  Count = c(11820, 9650, 8240, 7060, 6480, 10230, 1780)
)
annotation_flower(
  annotation_summary, total = 14200, style = style,
  output_file = "docs/images/transcriptome/Annotation_flower_sanshu.pdf"
)

function_summary <- data.frame(
  Category = c("Signal transduction", "Carbohydrate metabolism", "Translation",
               "Protein turnover", "Lipid metabolism", "Cell cycle",
               "Secondary metabolism", "Defense mechanisms"),
  Count = c(1260, 1120, 930, 810, 650, 520, 470, 390)
)
annotation_bar(
  function_summary, top_n = 8, title = "Functional classification",
  style = style, output_file = "docs/images/transcriptome/Annotation_bar_sanshu.pdf"
)

species_summary <- data.frame(
  Category = c("Species alpha", "Species beta", "Species gamma", "Species delta",
               "Species epsilon", "Species zeta", "Species eta", "Species theta",
               "Species iota", "Species kappa", "Unclassified"),
  Count = c(3820, 2460, 1740, 1290, 930, 710, 540, 420, 310, 240, 870)
)
annotation_pie(
  species_summary, top_n = 8, title = "NR species distribution",
  style = style, output_file = "docs/images/transcriptome/Annotation_pie_sanshu.pdf"
)
```

结果：

- [注释花瓣图](../images/transcriptome/Annotation_flower_sanshu.pdf)
- [功能分类柱状图](../images/transcriptome/Annotation_bar_sanshu.pdf)
- [物种组成饼图](../images/transcriptome/Annotation_pie_sanshu.pdf)

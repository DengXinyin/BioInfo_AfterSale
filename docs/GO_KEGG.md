# GO/KEGG 富集气泡图

本教程使用虚构的富集结果，复现参考图的布局：横轴为 `GeneRatio` 或 `RichFactor`，气泡大小为 `Gene count`，颜色图例在使用 P 值时为 `-log10(Pvalue)` 或 `-log10(Padj)`。所有连续图例只显示 3 个刻度文字，示例只输出 PDF。

## 准备数据

```r
library(BioInfoAfterSale)
set.seed(20260904)
enrichment_table <- data.frame(
  Description = c("Epidermis development", "Skin development", "Inflammatory response", "Cell differentiation", "Keratinocyte differentiation", "Keratinization", "Cornification", "Leukocyte chemotaxis", "Neutrophil chemotaxis", "Peptide cross-linking"),
  pvalue = c(0.00008, 0.0002, 0.0007, 0.0015, 0.003, 0.006, 0.009, 0.014, 0.021, 0.032),
  p.adjust = c(0.0010, 0.0018, 0.0045, 0.0080, 0.012, 0.019, 0.026, 0.034, 0.041, 0.048),
  GeneRatio = c(.065, .060, .050, .048, .038, .030, .026, .023, .020, .015),
  RichFactor = c(.64, .58, .51, .47, .42, .36, .31, .27, .22, .18),
  Count = c(60, 50, 45, 40, 35, 30, 25, 20, 15, 10), stringsAsFactors = FALSE
)
style <- choose_plot_style(font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 9, figure_height = 6.5, title = list(size = 18, bold = TRUE),
  axis_title = list(size = 15, bold = TRUE), axis_text = list(size = 12, bold = TRUE),
  legend_title = list(size = 14, bold = TRUE), legend_text = list(size = 12),
  legend = list(position = "right", frame = FALSE))
dir.create("docs/images", recursive = TRUE, showWarnings = FALSE)
```

## GeneRatio–Gene count–`-log10(Pvalue)`

```r
p_value <- GO_KEGG_plot(enrichment_table, plot_type = "dotplot", filter_by = "p.adjust", cutoff = .05,
  show_category = 10, x = "GeneRatio", color = "pvalue", x_label = "GeneRatio",
  color_transform = "neg_log10", color_label = expression(-log[10](Pvalue)),
  color_breaks = c(1.5, 3, 4.0), size = "Count", size_breaks = c(10, 30, 60), size_label = "Gene count",
  color_palette = c("#496B8E", "#A7A1B7", "#D9565B"), title = "GO enrichment", style = style,
  output_file = "docs/images/GO_KEGG_pvalue.pdf")
p_value
```

结果：[GO_KEGG_pvalue.pdf](images/GO_KEGG_pvalue.pdf)

颜色值在绘图前计算为 `-log10(pvalue)`；颜色图例和气泡大小图例均只保留 3 个文字刻度。

## GeneRatio–Gene count–`-log10(Padj)`

```r
p_padj <- GO_KEGG_plot(enrichment_table, plot_type = "dotplot", filter_by = "p.adjust", cutoff = .05,
  show_category = 10, x = "GeneRatio", color = "p.adjust", x_label = "GeneRatio",
  color_transform = "neg_log10", color_label = expression(-log[10](Padj)), color_breaks = c(1.4, 2, 3),
  size = "Count", size_breaks = c(10, 30, 60), size_label = "Gene count",
  color_palette = c("#496B8E", "#A7A1B7", "#D9565B"), title = "GO enrichment with adjusted P value", style = style,
  output_file = "docs/images/GO_KEGG_padj.pdf")
p_padj
```

结果：[GO_KEGG_padj.pdf](images/GO_KEGG_padj.pdf)

## Rich factor–Gene count–`-log10(Padj)`

```r
p_rich <- GO_KEGG_plot(enrichment_table, plot_type = "dotplot", filter_by = "p.adjust", cutoff = .05,
  show_category = 10, x = "RichFactor", color = "p.adjust", x_label = "Rich factor",
  color_transform = "neg_log10", color_label = expression(-log[10](Padj)), color_breaks = c(1.4, 2, 3),
  size = "Count", size_breaks = c(10, 30, 60), size_label = "Gene count",
  color_palette = c("#496B8E", "#A7A1B7", "#D9565B"), title = "GO enrichment by Rich factor", style = style,
  output_file = "docs/images/GO_KEGG_richfactor.pdf")
p_rich
```

结果：[GO_KEGG_richfactor.pdf](images/GO_KEGG_richfactor.pdf)

`GeneRatio`/`RichFactor` 作为横轴时不作为图例；如果把它们映射到 `color`，也可通过 `color_breaks` 限制为 3 个图例文字。`size_breaks` 应传入实际 `Count` 范围内的 3 个数值。

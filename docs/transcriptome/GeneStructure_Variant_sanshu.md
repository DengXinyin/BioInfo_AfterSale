# 转录组基因结构与变异统计图

本教程演示可变剪接事件、差异外显子使用、SNP/INDEL 数量和变异功能区域分布。柱状图与饼图复用通用注释汇总函数，DEU 使用独立的表达曲线加基因模型。全部数据均为虚构数据。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 10, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 13),
  axis_text = list(size = 10), legend = list(position = "right", frame = FALSE)
)

as_events <- expand.grid(
  Category = c("Treatment_vs_Control", "Recovery_vs_Control"),
  Event = c("SE", "A5", "A3", "MX", "RI"), stringsAsFactors = FALSE
)
as_events$Count <- c(218, 164, 106, 84, 63, 176, 139, 92, 71, 54)
annotation_bar(
  as_events, category = "Category", value = "Count", fill = "Event",
  horizontal = FALSE, position = "dodge", title = "Alternative splicing events",
  style = style, output_file = "docs/images/transcriptome/GeneStructure_AS_sanshu.pdf"
)

deu <- expand.grid(
  Exon = paste0("E", 1:8), Group = c("Control", "Treatment"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
deu$Expression <- c(42, 55, 61, 48, 73, 66, 57, 50,
                    45, 58, 96, 51, 78, 70, 60, 53)
deu$Padj <- rep(c(0.42, 0.31, 0.001, 0.27, 0.08, 0.11, 0.38, 0.49), 2)
deu_exon_expression(
  deu, colors = c(Control = "#4472C4", Treatment = "#D62728"),
  style = style, output_file = "docs/images/transcriptome/GeneStructure_DEU_sanshu.pdf",
  width = 10, height = 7
)

variant_counts <- expand.grid(
  Category = paste0("RNA-", LETTERS[1:6]), Type = c("SNPs", "INDELs"),
  stringsAsFactors = FALSE
)
variant_counts$Count <- c(42800, 39100, 44600, 41200, 43700, 40500,
                          8200, 7600, 8500, 7900, 8300, 7700)
annotation_bar(
  variant_counts, category = "Category", value = "Count", fill = "Type",
  horizontal = FALSE, position = "dodge",
  colors = c(SNPs = "#4E7CA1", INDELs = "#D67E56"),
  title = "Variant counts", style = style,
  output_file = "docs/images/transcriptome/Variant_counts_sanshu.pdf"
)

variant_effect <- data.frame(
  Category = c("EXON", "INTRON", "INTERGENIC", "UPSTREAM", "DOWNSTREAM",
               "UTR_3_PRIME", "UTR_5_PRIME", "SPLICE_SITE_REGION"),
  Count = c(12800, 17400, 9600, 6300, 5900, 2100, 1800, 920)
)
annotation_pie(
  variant_effect, top_n = 8, title = "Variant effect distribution",
  style = style, output_file = "docs/images/transcriptome/Variant_effect_sanshu.pdf"
)
```

结果：

- [可变剪接事件统计](../images/transcriptome/GeneStructure_AS_sanshu.pdf)
- [差异外显子使用](../images/transcriptome/GeneStructure_DEU_sanshu.pdf)
- [SNP/INDEL 数量](../images/transcriptome/Variant_counts_sanshu.pdf)
- [变异功能区域分布](../images/transcriptome/Variant_effect_sanshu.pdf)

# 转录组富集分析扩展图

本教程补充上下调蝴蝶图、多轨富集圈图和 GO 有向无环图；常规条形图、气泡图与 GSEA 图继续使用 `GO_KEGG_plot()` 和 `plot_gsea()`。示例术语、基因数量、P 值及节点关系均为虚构数据。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 10, figure_height = 8,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 13),
  axis_text = list(size = 9), legend_title = list(size = 11),
  legend_text = list(size = 10), legend = list(position = "right", frame = FALSE)
)

terms <- data.frame(
  Term = paste0("GO:", sprintf("%07d", 2000001:2000012), ": process ", LETTERS[1:12]),
  Up = c(82, 76, 64, 58, 49, 44, 39, 35, 31, 27, 24, 20),
  Down = c(61, 69, 52, 47, 55, 40, 36, 29, 25, 30, 18, 16)
)
enrichment_butterfly(
  terms, style = style,
  output_file = "docs/images/transcriptome/Enrichment_butterfly_sanshu.pdf",
  width = 11, height = 8
)

circle_data <- data.frame(
  ID = paste0("GO:", sprintf("%07d", 3000001:3000015)),
  Category = rep(c("BP", "CC", "MF"), each = 5),
  Count = c(180, 145, 126, 108, 92, 160, 132, 115, 96, 81, 138, 119, 102, 88, 73),
  Up = c(72, 64, 55, 47, 39, 61, 54, 46, 41, 32, 58, 49, 43, 35, 29),
  Down = c(108, 81, 71, 61, 53, 99, 78, 69, 55, 49, 80, 70, 59, 53, 44),
  RichFactor = seq(0.52, 0.18, length.out = 15),
  Padj = 10^seq(-8, -2, length.out = 15)
)
enrichment_circle(
  circle_data,
  category_colors = c(BP = "#F7CC13", CC = "#954572", MF = "#0796E0"),
  style = style, output_file = "docs/images/transcriptome/Enrichment_circle_sanshu.pdf",
  width = 11, height = 11
)

dag_nodes <- data.frame(
  ID = paste0("GO:", sprintf("%07d", 4000001:4000007)),
  Term = c("root process", "response pathway", "metabolic pathway",
           "cellular response", "redox process", "transport process", "target process"),
  Pvalue = c(0.04, 0.012, 0.018, 0.004, 0.0012, 0.007, 0.0002),
  GeneRatio = c("18/500", "16/500", "15/500", "13/500", "11/500", "10/500", "8/500")
)
dag_edges <- data.frame(
  Parent = dag_nodes$ID[c(1, 1, 2, 2, 3, 4, 5)],
  Child = dag_nodes$ID[c(2, 3, 4, 5, 6, 7, 7)]
)
go_dag(
  dag_nodes, dag_edges, style = style,
  output_file = "docs/images/transcriptome/Enrichment_GO_DAG_sanshu.pdf",
  width = 12, height = 8
)
```

结果：

- [上下调蝴蝶图](../images/transcriptome/Enrichment_butterfly_sanshu.pdf)
- [多轨富集圈图](../images/transcriptome/Enrichment_circle_sanshu.pdf)
- [GO DAG](../images/transcriptome/Enrichment_GO_DAG_sanshu.pdf)

KEGG 原生 pathway map 需要对应物种、版本一致的 KGML 和底图，不能只根据富集汇总表可靠重建；输入缺失时应明确停止，不生成说明页冒充结果图。

# 转录组差异表达汇总图

本教程覆盖差异基因数量、比较组花瓣图、UpSet 交集、MA 图、火山图和表达趋势。聚类热图继续使用通用 `Heatmap_plot()`。所有示例基因 ID 和统计量均为模拟数据。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 13),
  axis_text = list(size = 10), legend = list(position = "right", frame = FALSE)
)

deg_counts <- expand.grid(
  Comparison = c("Early_vs_Control", "Late_vs_Control", "Late_vs_Early"),
  Direction = c("Up", "Down"), stringsAsFactors = FALSE
)
deg_counts$Count <- c(620, 910, 480, 430, 760, 350)
deg_count_bar(
  deg_counts, style = style,
  output_file = "docs/images/transcriptome/DEG_counts_sanshu.pdf"
)

set.seed(2043)
universe <- paste0("Gene", sprintf("%04d", 1:1800))
sets <- list(
  Early_vs_Control = sample(universe, 620),
  Late_vs_Control = sample(universe, 910),
  Late_vs_Early = sample(universe, 480),
  Recovery_vs_Late = sample(universe, 540)
)
deg_flower(
  sets, style = style,
  output_file = "docs/images/transcriptome/DEG_flower_sanshu.pdf"
)
deg_upset(
  sets, max_intersections = 16, style = style,
  output_file = "docs/images/transcriptome/DEG_upset_sanshu.pdf",
  width = 11, height = 7
)

n_gene <- 2500
deg_table <- data.frame(
  MeanExpression = 2^runif(n_gene, 0, 15),
  log2FC = rnorm(n_gene, 0, 1.7),
  pvalue = runif(n_gene), stringsAsFactors = FALSE
)
deg_table$Status <- ifelse(deg_table$pvalue < 0.05 & deg_table$log2FC > 1, "Up",
  ifelse(deg_table$pvalue < 0.05 & deg_table$log2FC < -1, "Down", "NoSig"))
deg_ma(
  deg_table, style = style,
  output_file = "docs/images/transcriptome/DEG_MA_sanshu.pdf"
)
plot_volcano(
  transform(deg_table, log2FC = log2FC), log2fc = "log2FC", pvalue = "pvalue",
  status = "Status", status_levels = c("Up", "NoSig", "Down"),
  status_colors = c(Up = "#D5695D", NoSig = "#BDBDBD", Down = "#65A479"),
  style = style, output_file = "docs/images/transcriptome/DEG_volcano_sanshu.pdf"
)

trend <- expand.grid(
  Gene = paste0("TrendGene", 1:90), Sample = paste0("Time", 0:5),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
trend$Cluster <- rep(rep(c("Cluster 1", "Cluster 2", "Cluster 3"), each = 30), times = 6)
time <- rep(0:5, each = 90)
trend$Zscore <- rnorm(nrow(trend), 0, 0.3) + ifelse(
  trend$Cluster == "Cluster 1", time / 2.5,
  ifelse(trend$Cluster == "Cluster 2", -time / 2.5, sin(time / 1.5))
)
expression_trend(
  trend, sample_order = paste0("Time", 0:5), style = style,
  output_file = "docs/images/transcriptome/DEG_trend_sanshu.pdf",
  width = 11, height = 7
)
```

结果：

- [差异基因数量](../images/transcriptome/DEG_counts_sanshu.pdf)
- [比较组花瓣图](../images/transcriptome/DEG_flower_sanshu.pdf)
- [UpSet 交集图](../images/transcriptome/DEG_upset_sanshu.pdf)
- [MA 图](../images/transcriptome/DEG_MA_sanshu.pdf)
- [火山图](../images/transcriptome/DEG_volcano_sanshu.pdf)
- [表达趋势图](../images/transcriptome/DEG_trend_sanshu.pdf)

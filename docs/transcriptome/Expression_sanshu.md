# 转录组表达量汇总图

表达量分析复用包内通用箱线图、分布图和 PCA，仅新增带数值与分组标签的样本相关性热图。以下均为模拟表达量；矩阵行是基因、列是样本，样本信息通过同名标识严格匹配。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 13),
  axis_text = list(size = 10), legend = list(position = "right", frame = FALSE),
  group_palette = c("#4472C4", "#ED7D31", "#70AD47")
)

set.seed(2042)
samples <- paste0(rep(c("Control", "Early", "Late"), each = 4), "-", 1:4)
groups <- setNames(rep(c("Control", "Early", "Late"), each = 4), samples)
expression <- matrix(rlnorm(900 * length(samples), 2.1, 0.8), nrow = 900,
                     dimnames = list(paste0("Gene", seq_len(900)), samples))
expression[1:120, groups == "Late"] <- expression[1:120, groups == "Late"] * 2.2
long_expression <- data.frame(
  Sample = rep(samples, each = nrow(expression)),
  Group = rep(groups[samples], each = nrow(expression)),
  log2FPKM = as.vector(log2(expression + 1))
)

plot_violin_box(
  long_expression, value_column = "log2FPKM", group_column = "Sample",
  fill_column = "Group", plot_type = "box", show_jitter = FALSE,
  style = style, output_file = "docs/images/transcriptome/Expression_box_sanshu.pdf",
  figure_width = 10, figure_height = 6
)
plot_distribution(
  long_expression, value_column = "log2FPKM", group_column = "Sample",
  plot_type = "density", style = style,
  output_file = "docs/images/transcriptome/Expression_density_sanshu.pdf",
  figure_width = 10, figure_height = 7
)

correlation <- stats::cor(log2(expression + 1), method = "pearson")
sample_correlation_heatmap(
  correlation, group = groups, cluster = FALSE, digits = 2, style = style,
  output_file = "docs/images/transcriptome/Expression_correlation_sanshu.pdf",
  width = 9, height = 8
)

pca_data <- data.frame(
  Sample = samples, Group = unname(groups[samples]),
  t(log2(expression + 1)), check.names = FALSE
)
plot_pca(
  pca_data, expr_columns = rownames(expression), show_ellipse = TRUE,
  show_labels = TRUE, style = style,
  output_file = "docs/images/transcriptome/Expression_PCA_sanshu.pdf"
)
```

结果：

- [表达量箱线图](../images/transcriptome/Expression_box_sanshu.pdf)
- [表达量密度图](../images/transcriptome/Expression_density_sanshu.pdf)
- [样本相关性热图](../images/transcriptome/Expression_correlation_sanshu.pdf)
- [表达量 PCA](../images/transcriptome/Expression_PCA_sanshu.pdf)

# WGCNA 结果可视化补充

已有 `WGCNA_plot_soft_threshold()`、`Heatmap_plot()` 和 `ppi_network()` 分别覆盖软阈值、模块基因热图及网络图。本教程补充样本树、模块基因数量、模块—性状热图和 MM–GS 散点图。数据均为模拟结果，不代表真实生物学结论。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 9, figure_height = 7,
  title = list(size = 18, bold = TRUE), axis_title = list(size = 13),
  axis_text = list(size = 10), legend = list(position = "right", frame = FALSE)
)

set.seed(2045)
sample_expression <- matrix(
  rnorm(12 * 150), nrow = 12,
  dimnames = list(paste0("Sample", sprintf("%02d", 1:12)), paste0("Gene", 1:150))
)
wgcna_sample_tree(
  sample_expression, style = style,
  output_file = "docs/images/transcriptome/WGCNA_sample_tree_sanshu.pdf"
)

module_sizes <- data.frame(
  Module = c("turquoise", "blue", "brown", "yellow", "green", "red"),
  Count = c(1860, 1320, 940, 710, 530, 340)
)
wgcna_module_sizes(
  module_sizes, style = style,
  output_file = "docs/images/transcriptome/WGCNA_module_sizes_sanshu.pdf"
)

module_trait_cor <- matrix(
  c(0.72, -0.31, 0.18, -0.64, 0.55, 0.21,
    -0.48, 0.67, -0.12, 0.33, -0.58, 0.44),
  nrow = 6, dimnames = list(module_sizes$Module, c("Treatment", "Time"))
)
module_trait_p <- matrix(
  c(0.004, 0.18, 0.42, 0.011, 0.036, 0.35,
    0.07, 0.008, 0.61, 0.20, 0.028, 0.09),
  nrow = 6, dimnames = dimnames(module_trait_cor)
)
wgcna_module_trait_heatmap(
  module_trait_cor, module_trait_p, style = style,
  output_file = "docs/images/transcriptome/WGCNA_module_trait_sanshu.pdf"
)

mm_gs <- data.frame(
  ModuleMembership = pmin(1, pmax(0, rbeta(260, 4, 2))),
  Module = rep(c("turquoise", "blue"), each = 130)
)
mm_gs$GeneSignificance <- pmin(1, pmax(0,
  0.72 * mm_gs$ModuleMembership + rnorm(nrow(mm_gs), 0, 0.13)))
wgcna_mm_gs(
  mm_gs, module = "Module",
  colors = c(turquoise = "turquoise3", blue = "steelblue"),
  style = style, output_file = "docs/images/transcriptome/WGCNA_MM_GS_sanshu.pdf"
)
```

结果：

- [样本聚类树](../images/transcriptome/WGCNA_sample_tree_sanshu.pdf)
- [模块基因数量](../images/transcriptome/WGCNA_module_sizes_sanshu.pdf)
- [模块—性状热图](../images/transcriptome/WGCNA_module_trait_sanshu.pdf)
- [MM–GS 散点图](../images/transcriptome/WGCNA_MM_GS_sanshu.pdf)

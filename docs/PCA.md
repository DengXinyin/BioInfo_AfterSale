# PCA 可视化示例

使用 `plot_pca()` 绘制样本主成分分析图。函数兼容 `choose_plot_style()`，可统一设置字体、字号、图例、面板、画布尺寸和 DPI。

所有字体优先使用 `Times New Roman`。PCA 图例放在主图内部且不带方框，样本标签使用避让算法减少重叠。

## 示例数据与绘图

```r
library(BioInfoAfterSale)

after_sale_style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 400,
  figure_width = 9, figure_height = 7,
  title = list(size = 24, bold = TRUE),
  axis_title = list(size = 18, bold = TRUE),
  axis_text = list(size = 16, bold = TRUE),
  legend_title = list(size = 16, bold = TRUE),
  legend_text = list(size = 15, bold = TRUE),
  legend = list(position = "right", frame = FALSE),
  panel = list(border = TRUE, major_grid = TRUE, minor_grid = FALSE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)

set.seed(2026)
pca_data <- data.frame(
  Sample = paste0(rep(c("Control", "TA0", "TA4"), each = 5), "-", 1:5),
  PC1 = c(rnorm(5, 0, 12), rnorm(5, 38, 13), rnorm(5, -34, 15)),
  PC2 = c(rnorm(5, -5, 8), rnorm(5, 2, 9), rnorm(5, 15, 10)),
  Group = rep(c("Control", "TA0%", "TA4%"), each = 5)
)

p_pca <- plot_pca(
  data = pca_data, pc1 = "PC1", pc2 = "PC2", group = "Group",
  sample = "Sample", variance = c(77.87, 7.60), ellipse = TRUE,
  show_labels = TRUE, label_repel = TRUE,
  title = "Principal Component Analysis", legend_inside = c(0.98, 0.98),
  style = after_sale_style, output_file = "docs/images/PCA.pdf"
)
p_pca
```

不显示样本标签时设置 `show_labels = FALSE`。

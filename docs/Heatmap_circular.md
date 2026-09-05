# 环形聚类热图教程

本示例使用**虚构**的基因表达矩阵（50 个基因 × 9 个样本），演示用
`plot_circular_heatmap()` 绘制环形（环状）聚类热图。参考的绘图思路来自
"如何使用 DeepSeek 绘制好看的环形热图"教程：对表达矩阵按行做 Z-score
归一化，用 `circlize::circos.heatmap()` 把基因画在圆周方向、样本画在径向，
基因按表达模式聚类后树状图显示在圆环内侧，基因名标注在外侧，并添加
径向的样本标签与连续色阶图例。示例只输出 PDF。

```r
library(BioInfoAfterSale)
set.seed(20260903)

# 虚构数据：50 个基因（行）x 9 个样本（列），模拟 RNA-seq 表达量
gene_ids <- paste0("Gene_", sprintf("%02d", 1:50))
sample_ids <- paste0("Sample_", 1:9)
expression <- matrix(
  rnorm(length(gene_ids) * length(sample_ids), mean = 8, sd = 1.5),
  nrow = length(gene_ids),
  dimnames = list(gene_ids, sample_ids)
)

style <- choose_plot_style(
  font_family = "Times New Roman",
  theme = "classic",
  dpi = 300,
  figure_width = 10,
  figure_height = 10,
  title = list(size = 18, bold = TRUE),
  legend_title = list(size = 14),
  legend_text = list(size = 12)
)
dir.create("docs/images", showWarnings = FALSE)

plot_circular_heatmap(
  expression,
  scale = "row",          # 按基因（行）Z-score 归一化
  cluster = TRUE,         # 基因聚类，聚类树显示在圆环内侧
  dend_side = "inside",
  show_rownames = TRUE,   # 基因名标在外侧
  rownames_side = "outside",
  rownames_cex = 0.7,
  show_colnames = TRUE,   # 样本名沿径向标注
  colnames_cex = 0.7,
  start_degree = 45,      # 起始角度
  gap_degree = 45,        # 扇区间隔
  heatmap_colors = c("#A5CC26", "white", "#FF7BAC"),
  legend_title = "Exp",
  title = "Circular expression heatmap",
  style = style,
  output_file = "docs/images/Heatmap_circular.pdf"
)
```

结果：[Heatmap_circular.pdf](images/Heatmap_circular.pdf)

图中每个基因（行）占据圆环上的一个位置，颜色表示该基因在对应样本中的
Z-score：绿色为低表达、白色为中位、粉色为高表达。相邻基因颜色模式相近时
会聚到同一分支，圆环内侧的树状图展示基因聚类关系；圆环外侧为基因名，
径向刻度附近为样本名。数据为随机模拟，仅用于演示函数用法，不代表真实
生物学实验。

## 常用参数调整

- 换配色：改 `heatmap_colors = c("#2166AC", "#F7F7F7", "#B2182B")`（蓝白红），
  同时可配合 `zscore_breaks = c(-2, 0, 2)` 调整颜色断点。
- 聚类树放外侧、基因名放内侧：`dend_side = "outside"`，
  `rownames_side = "inside"`。
- 样本数很多时不标注基因名：`show_rownames = FALSE`；
  基因数很多时不标注样本：`show_colnames = FALSE`，避免文字重叠。
- 需要位图时改 `output_file` 后缀为 `.png` 并设置 `dpi`。

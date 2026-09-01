# PCA、火山图、Venn 和 GSEA 可视化示例

所有函数均兼容 `choose_plot_style()`，可统一设置字体、字号、图例、面板、画布尺寸和 DPI。

本页的模拟数据用于快速理解函数接口；文末的“对应项目结果”一节展示了如何读取本项目
`Input/` 中的 PCA、差异基因、Venn 汇总和 GSEA 结果。函数只负责绘图，不会修改原始表格。

## 公共样式

```r
library(BioInfoAfterSale)

after_sale_style <- choose_plot_style(
  font_family = "Times New Roman",
  theme = "bw",
  dpi = 400,
  figure_width = 9,
  figure_height = 7,
  title = list(size = 18, bold = TRUE),
  axis_title = list(size = 14, bold = TRUE),
  axis_text = list(size = 12, bold = TRUE),
  legend_title = list(size = 13, bold = TRUE),
  legend_text = list(size = 12, bold = TRUE),
  data_label = list(size = 11, bold = TRUE),
  legend = list(position = "right", frame = FALSE),
  panel = list(border = TRUE, major_grid = TRUE, minor_grid = FALSE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)
```

## PCA

```r
set.seed(2026)
pca_data <- data.frame(
  Sample = paste0(rep(c("Control", "TA0", "TA4"), each = 5), "-", 1:5),
  PC1 = c(rnorm(5, 0, 12), rnorm(5, 38, 13), rnorm(5, -34, 15)),
  PC2 = c(rnorm(5, -5, 8), rnorm(5, 2, 9), rnorm(5, 15, 10)),
  Group = rep(c("Control", "TA0%", "TA4%"), each = 5)
)

p_pca <- plot_pca(
  data = pca_data,
  pc1 = "PC1", pc2 = "PC2", group = "Group", sample = "Sample",
  variance = c(77.87, 7.60), ellipse = TRUE, show_labels = TRUE,
  title = "Principal Component Analysis",
  style = after_sale_style,
  output_file = "PCA.pdf"
)
p_pca
```

不显示样本标签时设置 `show_labels = FALSE`。

## 火山图

```r
set.seed(2026)
n_gene <- 5000
deg_data <- data.frame(
  GeneID = paste0("Gene", seq_len(n_gene)),
  log2FC = rnorm(n_gene, 0, 1.8),
  pvalue = pmin(runif(n_gene), 10^(-rexp(n_gene, rate = 0.8)))
)

p_volcano <- plot_volcano(
  data = deg_data,
  log2fc = "log2FC", pvalue = "pvalue",
  fc_cutoff = 1, p_cutoff = 0.05,
  status_colors = c("#E25659", "#D4D4D4", "#335372"),
  title = "TA4% vs Control", legend_title = NULL,
  style = after_sale_style,
  output_file = "Volcano_TA4_vs_Control.png"
)
p_volcano
```

已有差异状态列时：

```r
deg_data$level <- ifelse(
  deg_data$log2FC >= 1 & deg_data$pvalue < 0.05, "Up",
  ifelse(deg_data$log2FC <= -1 & deg_data$pvalue < 0.05,
         "Down", "Not significant")
)
plot_volcano(deg_data, status = "level", legend_title = NULL,
             style = after_sale_style)
```

## 三集合 Venn

```r
set.seed(2026)
gene_universe <- paste0("Gene", 1:5000)
venn_sets <- list(
  `TA4% vs TA0%` = sample(gene_universe, 1800),
  `TA4% vs Control` = sample(gene_universe, 2200),
  `TA0% vs Control` = sample(gene_universe, 1500)
)

p_venn <- plot_venn(
  sets = venn_sets,
  colors = c("#A3A4CA", "#E25659", "#335372"),
  title = "Differentially Expressed Gene Overlap",
  style = after_sale_style,
  output_file = "Venn_DEG.pdf", width = 8, height = 7
)
p_venn
```

## GSEA

`hits` 可为命中位置或与排序向量等长的逻辑向量。

```r
set.seed(2026)
n_rank <- 12000
ranked_metric <- sort(rnorm(n_rank), decreasing = TRUE)
hit_positions <- sort(sample(seq_len(n_rank), 260))
is_hit <- seq_len(n_rank) %in% hit_positions
hit_weight <- abs(ranked_metric) * is_hit
running_score <- cumsum(
  hit_weight / sum(hit_weight) - (!is_hit) / sum(!is_hit)
)

p_gsea <- plot_gsea(
  running_score = running_score,
  hits = hit_positions,
  ranked_metric = ranked_metric,
  title = "Peptide cross-linking",
  statistics = c(NES = 2.22, `P value` = 6.1e-09,
                 `Adjusted P` = 2.8e-06),
  colors = c("#E25659", "#335372"),
  style = after_sale_style,
  output_file = "GSEA_peptide_cross_linking.pdf",
  width = 9.2, height = 6.8
)
p_gsea
```

## 改变图例位置

```r
left_legend_style <- choose_plot_style(
  font_family = "Times New Roman",
  legend = list(position = "left", frame = FALSE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)

plot_volcano(deg_data, legend_title = NULL, style = left_legend_style)
```

## 输出格式

`output_file` 根据扩展名选择输出格式：

```r
plot_pca(pca_data, style = after_sale_style, output_file = "PCA.pdf")
plot_pca(pca_data, style = after_sale_style, output_file = "PCA.png")
plot_pca(pca_data, style = after_sale_style, output_file = "PCA.svg")
```

不提供 `output_file` 时，函数返回 ggplot 对象，便于继续添加 ggplot2 图层。

## 对应项目结果

如果在本仓库目录结构中运行以下代码，`input_dir` 默认指向项目的 `Input/` 目录。将
`input_dir` 改成自己的数据目录后，也可以复用这套代码。

```r
library(readxl)

input_dir <- "../Input"  # 从 BioInfo_AfterSale/ 目录运行
if (!dir.exists(input_dir)) {
  stop("请将 input_dir 修改为包含项目 Excel 文件的目录。")
}
```

### PCA：`all.PCA.xlsx`

该表第一列是样本名，最后一行是解释方差比例，不应作为样本绘图。这里从样本名中
提取 `Control`、`TA0` 和 `TA4` 作为分组，并使用前两列主成分。

```r
pca_project <- readxl::read_excel(file.path(input_dir, "all.PCA.xlsx"))
names(pca_project)[1] <- "Sample"
pca_project <- pca_project[pca_project$Sample != "Explained Variance Ratio", ]
pca_project$Group <- sub("-[^-]+$", "", pca_project$Sample)

plot_pca(
  pca_project,
  pc1 = "PC1", pc2 = "PC2", group = "Group", sample = "Sample",
  variance = c(77.87, 7.60),
  style = after_sale_style,
  output_file = "PCA_project.pdf"
)
```

### 火山图：差异基因表

`B_vs_A.DEGs.xlsx` 已包含 `log2FC`、`pvalue` 和 `level`。项目表中的状态值为
`up`、`down`、`nosig`，因此先映射为绘图函数默认的图例名称。

```r
deg_project <- readxl::read_excel(file.path(input_dir, "B_vs_A.DEGs.xlsx"))
deg_project$status <- c(up = "Up", down = "Down", nosig = "Not significant")[
  tolower(deg_project$level)
]

plot_volcano(
  deg_project,
  log2fc = "log2FC", pvalue = "pvalue", status = "status",
  status_levels = c("Up", "Not significant", "Down"),
  status_colors = c("#E25659", "#CFCFCF", "#335372"),
  title = "B vs A", style = after_sale_style,
  output_file = "Volcano_B_vs_A.png"
)
```

缺失的 `log2FC` 或 `pvalue` 会在图中产生缺失点；正式出图前可根据项目规则筛除
这些行，例如 `deg_project <- subset(deg_project, !is.na(log2FC) & !is.na(pvalue))`。

### Venn：`Summary_Venn.xlsx`

Venn 汇总表用 0/1 列表示基因是否属于某个集合。`plot_venn()` 接收的是三个向量，
所以需要先按标记列筛选 `GeneID`。

```r
venn_project <- readxl::read_excel(file.path(input_dir, "Summary_Venn.xlsx"))
venn_sets <- lapply(
  venn_project[c("C_vs_B", "C_vs_A", "B_vs_A")],
  function(flag) venn_project$GeneID[flag == 1]
)

plot_venn(
  venn_sets,
  title = "Differentially Expressed Gene Overlap",
  style = after_sale_style,
  output_file = "Venn_project.pdf", width = 8, height = 7
)
```

### GSEA：使用排序向量和命中位置

`C_vs_A.GSEA.GOBP.xlsx`、`GOCC.xlsx` 和 `KEGG.xlsx` 保存了 NES、P 值及命中排名，
但不包含完整的排序统计量。因此，绘制 GSEA 曲线时仍需使用运行 GSEA 时的
`ranked_metric` 和基因集命中位置；不能仅凭结果表准确重建曲线。

下面的代码展示如何用结果表中的第一条通路填写统计信息。`ranked_metric` 应替换为
实际 GSEA 输入的、按降序排列的全基因排序向量，`hit_positions` 应替换为该通路的
命中位置。

```r
gsea_table <- readxl::read_excel(file.path(input_dir, "C_vs_A.GSEA.GOBP.xlsx"))
term <- gsea_table[1, ]

# 用真实 GSEA 输入替换这两个示例对象。
ranked_metric <- sort(rnorm(12000), decreasing = TRUE)
hit_positions <- sort(sample(seq_along(ranked_metric), 260))
is_hit <- seq_along(ranked_metric) %in% hit_positions
running_score <- cumsum(
  abs(ranked_metric) * is_hit / sum(abs(ranked_metric) * is_hit) -
    (!is_hit) / sum(!is_hit)
)

plot_gsea(
  running_score = running_score, hits = hit_positions,
  ranked_metric = ranked_metric, title = term$Description,
  statistics = c(NES = term$NES, `P value` = term$pvalue,
                 `Adjusted P` = term$p.adjust),
  style = after_sale_style,
  output_file = "GSEA_project.pdf", width = 9.2, height = 6.8
)
```

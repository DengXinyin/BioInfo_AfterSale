# 用 `enrichment_bar()` 绘制富集分析条形图

本教程中的所有字体优先使用 `Times New Roman`。

本教程演示如何用 `BioInfoAfterSale::enrichment_bar()` 绘制带分组色带、Gene
count 气泡、通路名称和基因标签的富集分析条形图。示例数据直接写在文档中，
复制代码即可运行，不需要额外准备 Excel 文件。

> 示例数据只用于演示绘图接口，不代表真实 KEGG 或 GO 富集结果。

## 准备工作

首次使用时从 GitHub 安装 R 包：

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("DengXinyin/BioInfo_AfterSale")
```

载入包并创建结果目录：

```r
library(BioInfoAfterSale)

dir.create("enrichment_bar_results", recursive = TRUE, showWarnings = FALSE)
```

## 第 1 步：准备示例富集结果

`enrichment_bar()` 至少需要 5 列：

- `group`：分组名称，例如细胞类型、比较组或数据库类别；
- `Description`：富集条目或通路名称；
- `pvalue`：P 值，绘图时转换为 `-log10(pvalue)`；
- `count`：命中基因数，用气泡大小表示；
- `geneID`：命中基因，多个基因可用 `/` 分隔。

```r
enrichment_table <- data.frame(
  group = c(
    rep("Melanocyte", 6),
    rep("Keratinocyte", 6)
  ),
  Description = c(
    "Hypertrophic cardiomyopathy",
    "Dilated cardiomyopathy",
    "Leukocyte transendothelial migration",
    "cGMP-PKG signaling pathway",
    "Human papillomavirus infection",
    "Adrenergic signaling in cardiomyocytes",
    "Tuberculosis",
    "TGF-beta signaling pathway",
    "NOD-like receptor signaling pathway",
    "Malaria",
    "Hepatitis B",
    "Hepatitis C"
  ),
  pvalue = c(
    1.36e-9, 1.58e-7, 4.25e-6, 4.72e-5, 5.29e-4, 1.34e-3,
    1.00e-6, 4.72e-7, 1.00e-5, 1.71e-5, 1.94e-4, 2.00e-3
  ),
  count = c(15, 20, 21, 30, 15, 16, 40, 121, 16, 18, 30, 9),
  geneID = c(
    "Lama2/Tgfb2/Ttn/Myh7/Slc8a1",
    "Lama2/Tgfb2/Ttn/Myh7/Slc8a1",
    "Ncf1/Thy1/Cldn18/Cxcr4/Ptk2b",
    "Nppb/Atp2b2/Adrb3/Mylk3/Myh7/Slc8a1",
    "Atp6v1b1/Lama2/Thbs2/Ppp2r2b/Tnn/Wnt8b",
    "Ppp2r2b/Camk2a/Atp2b2/Myh7/Slc8a1",
    "Mrc2/Atp6v0a4/Tgfb2/Tlr1",
    "Ltbp1/Tgfb2/Inhba",
    "Trpm2/Irf7/Pstpip1/Oas3",
    "Thbs2/Tgfb2",
    "Irf7/Tgfb2/Ptk2b",
    "Ppp2r2b/Irf7/Oas3"
  ),
  stringsAsFactors = FALSE
)
```

## 第 2 步：创建统一绘图样式

`enrichment_bar()` 使用 `choose_plot_style()` 管理颜色、字体、字号、图例和
面板主题。这里使用 `classic` 主题，并给两个分组指定颜色。

```r
bar_style <- choose_plot_style(
  font_family = "Times New Roman",
  theme = "classic",
  dpi = 300,
  figure_width = 8,
  figure_height = 8,
  axis_title = list(size = 19),
  axis_text = list(size = 12),
  data_label = list(size = 10),
  facet_label = list(size = 14, bold = TRUE),
  legend_title = list(size = 15, bold = TRUE),
  legend_text = list(size = 12),
  legend = list(position = "right", frame = FALSE),
  panel = list(
    border = FALSE,
    major_grid = FALSE,
    minor_grid = FALSE
  ),
  group_palette = c(
    Melanocyte = "#F5BA68",
    Keratinocyte = "#8E9985"
  )
)
```

## 第 3 步：绘制富集条形图

```r
p <- enrichment_bar(
  data = enrichment_table,
  style = bar_style,
  output_file = "docs/images/enrichment_bar.png"
)

p
```

运行后会得到：

```r
file.exists("enrichment_bar_results/enrichment_bar.png")
file.info("enrichment_bar_results/enrichment_bar.png")$size
```

图中横轴为 `-log10(pvalue)`，条形越长表示 P 值越小；左侧圆圈大小和圆圈中的
数字表示 `count`；条形和左侧分组色带使用 `group_palette` 中的颜色。

## 第 4 步：隐藏基因标签

如果通路较多或基因名称太长，可以设置 `gene_scale = 0` 隐藏基因标签：

```r
enrichment_bar(
  data = enrichment_table,
  style = bar_style,
  gene_scale = 0,
  output_file = "docs/images/enrichment_bar_without_genes.png"
)
```

## 第 5 步：调整图例或文字

图例位置由 `choose_plot_style()` 控制。例如把图例放到底部，并适当减小文字：

```r
bottom_style <- choose_plot_style(
  font_family = "Times New Roman",
  theme = "classic",
  dpi = 300,
  figure_width = 8,
  figure_height = 8,
  axis_title = list(size = 18),
  axis_text = list(size = 11),
  data_label = list(size = 9),
  facet_label = list(size = 13, bold = TRUE),
  legend_title = list(size = 13, bold = TRUE),
  legend_text = list(size = 11),
  legend = list(position = "bottom", frame = FALSE, box = "horizontal"),
  panel = list(border = FALSE, major_grid = FALSE, minor_grid = FALSE),
  group_palette = c(
    Melanocyte = "#F5BA68",
    Keratinocyte = "#8E9985"
  )
)

enrichment_bar(
  data = enrichment_table,
  style = bottom_style,
  output_file = "docs/images/enrichment_bar_bottom_legend.png"
)
```

## 使用自己的数据

实际分析时，只要把自己的富集结果整理为同样的 5 列即可：

```r
my_table <- readxl::read_excel("my_enrichment_result.xlsx")

enrichment_bar(
  data = my_table,
  style = bar_style,
  output_file = "docs/images/my_enrichment_bar.png"
)
```

`pvalue` 必须是大于 0 的数值；缺失值、非数值和 `pvalue <= 0` 的记录会被过滤。

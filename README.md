# BioInfoAfterSale

`BioInfoAfterSale` 是面向生物信息售后分析的内部 R 包。第一版提供 GO/KEGG
过度富集分析及可视化，统计计算基于 `clusterProfiler`，基础绘图基于
`enrichplot`，并使用 `ggplot2` 统一常见的字体、字号、标题和图例设置。

## 主要函数

- `GO_KEGG_analyse()`：只做 GO 或 KEGG 富集分析；
- `GO_KEGG_plot()`：只对已有的 `enrichResult`/`gseaResult` 绘图；
- `GO_KEGG()`：统一入口，自动串联前两个函数。
- `choose_plot_style()`：创建可供所有 ggplot2 图形复用的统一样式。
- `Heatmap_plot()`：使用 ComplexHeatmap 绘制带 Z-score 和分组图例的热图。

## 安装

```r
install.packages("remotes")
remotes::install_github("DengXinyin/BioInfo_AfterSale")
```

Docker 镜像中应提前安装 `clusterProfiler`、`enrichplot`、`ggplot2`，以及实际
使用的物种注释包。

## 一键完成 GO 富集和绘图

```r
library(BioInfoAfterSale)

go <- GO_KEGG(
  gene = gene_ids,
  analysis = "GO",
  run_args = list(
    species = "human",
    key_type = "ENTREZID",
    ont = "ALL",
    pvalue_cutoff = 0.05
  ),
  plot_args = list(
    filter_by = "pvalue",
    cutoff = 0.05,
    color = "pvalue",
    show_category = 20,
    title = "GO Enrichment",
    font_family = "Arial",
    base_size = 16
  )
)

go$plot
go$table
```

## 统一可视化样式

```r
style <- choose_plot_style(
  font_family = "Arial",
  theme = "bw",
  dpi = 300,
  figure_width = 12,
  figure_height = 10,
  title = list(size = 22, bold = TRUE, align = "center"),
  axis_title = list(size = 18),
  axis_text = list(size = 16),
  legend_title = list(size = 18),
  legend_text = list(size = 15),
  facet_label = list(size = 17, bold = TRUE),
  legend = list(
    show = TRUE, position = "right", frame = FALSE,
    box = "vertical", box_just = "center"
  ),
  panel = list(
    border = TRUE,
    border_color = "black",
    border_width = 0.8,
    major_grid = TRUE,
    minor_grid = TRUE,
    major_grid_color = "#BDBDBD",
    minor_grid_color = "#E1E1E1"
  )
)

p <- GO_KEGG_plot(
  existing_enrich_result,
  filter_by = "pvalue",
  cutoff = 0.05,
  style = style
)
```

`GO_KEGG_plot()` 可以直接按英寸保存指定画布大小：

```r
p <- GO_KEGG_plot(
  selected_pathways,
  filter_by = NULL,
  x = "pvalue",
  x_transform = "neg_log10",
  x_label = expression(-log[10](Pvalue)),
  y_label = NULL,
  x_limits = c(1, 6),
  x_breaks = 1:6,
  x_expand = c(0.02, 0),
  color = "RichFactor",
  color_label = "Rich Factor",
  size = "Count",
  size_label = "Count",
  size_breaks = c(4, 7, 10),
  size_range = c(4, 12),
  point_alpha = 0.9,
  legend_order = c("color", "size"),
  output_file = "GO_enrichment.pdf",
  figure_width = 15,
  figure_height = 7,
  dpi = 300
)
```

不提供 `figure_width`、`figure_height` 或 `dpi` 时，使用
`choose_plot_style()` 样式对象中保存的默认值。

样式对象还保存了 `dpi`、默认画布宽高和 `group_palette`，后续其他售后绘图
函数可以直接复用同一套接口。

普通富集结果表也可以直接绘图。例如横轴使用 `-log10(pvalue)`、颜色使用
Rich Factor：

```r
p <- GO_KEGG_plot(
  selected_pathways,
  filter_by = NULL,
  x = "pvalue",
  x_transform = "neg_log10",
  color = "RichFactor",
  label = "Description",
  size = "Count",
  style = style
)
```

`pvalue_cutoff` 是富集计算阶段传给 `clusterProfiler` 的参数；`filter_by` 和
`cutoff` 是绘图前采用的筛选口径。因此可明确实现“按照原始 P 值而不是校正后
P 值筛选和绘图”。

## ComplexHeatmap 热图

```r
heatmap <- Heatmap_plot(
  expression_matrix,
  group = c(C1 = "Control", C2 = "Control", T1 = "Treatment", T2 = "Treatment"),
  group_colors = c(Control = "#7F8C8D", Treatment = "#8D6E63"),
  scale = "row",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_italic = TRUE,
  column_names_rot = 45,
  title = "Selected genes",
  title_position = "top",
  title_font_family = "Times New Roman",
  title_font_size = 18,
  title_font_face = "bold",
  legend_side = "right",
  output_file = "selected_genes_heatmap.pdf",
  figure_width = 8,
  figure_height = 7
)
```

默认把 Z-score 和 Group 图例放在主图右侧，Z-score 在上、Group 在下；
`column_names_rot = 0` 为水平，`90` 为竖直，也可传入任意角度。基因行名默认
使用斜体。`show_zscore_legend` 和 `show_group_legend` 可分别控制两类图例。

## 使用本地 KEGG 数据

离线 Docker 环境推荐提供固定版本的 `TERM2GENE` 和 `TERM2NAME`：

```r
kegg <- GO_KEGG(
  gene = entrez_ids,
  analysis = "KEGG",
  run_args = list(
    universe = background_entrez_ids,
    term2gene = kegg_term2gene,
    term2name = kegg_term2name
  ),
  plot_args = list(
    filter_by = "p.adjust",
    cutoff = 0.05,
    show_category = 15,
    title = "KEGG Enrichment"
  )
)
```

没有提供本地 `term2gene` 时，KEGG 会调用 `clusterProfiler::enrichKEGG()`；
这种方式可能需要访问 KEGG，因此不适合作为完全离线流程的默认方案。

## 已有富集结果时只重绘

```r
redrawn <- GO_KEGG(
  result = existing_enrich_result,
  plot_args = list(
    plot_type = "dotplot",
    filter_by = "pvalue",
    cutoff = 0.05,
    font_family = "Times New Roman",
    base_size = 18
  )
)
```

`GO_KEGG_plot()` 返回普通 ggplot 对象，仍可继续追加 `theme()`：

```r
p <- GO_KEGG_plot(existing_enrich_result, filter_by = "pvalue")
p + ggplot2::theme(legend.position = "bottom")
```

## 当前范围

版本 `0.2.0` 先覆盖 GO/KEGG ORA 的计算，以及 `dotplot` 和 `barplot`。指定
通路单独标红、GO 三分面、圈图、结果 Excel 导出等需求将在后续版本逐步加入。

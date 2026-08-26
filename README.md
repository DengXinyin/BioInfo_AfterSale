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
- `LEFSE_run()`：从 OTU 丰度、分类学表和分组信息计算 LEfSe marker。
- `LEFSE_LDA_plot()`：从 marker 表绘制 LDA score 柱状图。
- `LEFSE_cladogram_plot()`：绘制 LEfSe cladogram；支持两列分类学标签、单行分组图例和统一字体。
- `WGCNA_run()`：串联表达矩阵质控、soft power、模块、ME、kME和可选性状/网络分析。
- `WGCNA_prepare_expression()`：整理表达矩阵、变换、零值过滤和`goodSamplesGenes()`质控。
- `WGCNA_select_power()`：支持scale-free诊断优先或复刻镜像样本数规则。
- `WGCNA_module_sample()`：无外部表型时输出模块ME与样本/时间分组关系。

## WGCNA

WGCNA模块位于`R/WGCNA.R`，以
`192.168.30.202:23099/wgcna/wgcna:v1.6.7`中的`/script/wgcna2.r`为基准。
镜像的核心默认参数均可复现：`unsigned`网络/TOM、`deepSplit = 2`、
`minModuleSize = 30`、`mergeCutHeight = 0.25`、
`pamRespectsDendro = FALSE`、`maxBlockSize = 1000`及样本数soft power规则。

完整的实际数据教程见[`WGCNA_example_tutorial.md`](WGCNA_example_tutorial.md)，
可执行脚本见`inst/examples/WGCNA_example_tutorial.R`。

### R包依赖

WGCNA模块的非base直接依赖为`WGCNA`；其自身会使用`dynamicTreeCut`和
`fastcluster`。绘图和文件输出使用R自带的`stats`、`graphics`、`grDevices`、
`utils`和`tools`。每个函数体开头均以`# Required R packages:`注明对应依赖，
并使用`WGCNA::函数`进行调用。环境构建示例：

```r
# 在联网的R环境中执行一次；正式分析函数不会在运行过程中自动安装包
install.packages("WGCNA")
```

本包额外把流程拆成独立函数，并修正了两个不适合直接固化到通用包的行为：

- soft power默认优先选择达到scale-free阈值且斜率为负的首个候选；无候选时才回退镜像规则。需要严格复刻镜像时使用`strategy = "image_rule"`；
- 无表型时不构造“单个样本=1、其他样本=0”的伪性状，不输出缺乏统计解释的模块—单样本相关P值，而是输出每个模块在每个样本中的ME及可选分组汇总。

完整入口：

```r
expression <- read.delim("selected_gene_expression_fpkm.tsv", row.names = 1,
                         check.names = FALSE)
sample_info <- read.delim("sampleinfo.tsv", check.names = FALSE)

wgcna <- WGCNA_run(
  expression,
  sample_info = sample_info,
  time_col = "time_group",
  prepare_args = list(
    transform = "log2p1",
    keep_all_genes = TRUE,
    remove_bad = FALSE
  ),
  power_args = list(
    network_type = "unsigned",
    strategy = "scale_free",
    fit_cutoff = 0.9
  ),
  module_args = list(threads = 8),
  output_dir = "WGCNA_Output"
)
```

严格复刻v1.6.7镜像的power选择时：

```r
wgcna <- WGCNA_run(
  expression,
  sample_info = sample_info,
  prepare_args = list(transform = "none", keep_all_genes = TRUE),
  power_args = list(strategy = "image_rule", network_type = "unsigned"),
  module_args = list(threads = 72)
)
```

分步调用及主要返回内容：

```r
prepared <- WGCNA_prepare_expression(expression, transform = "log2p1")
power <- WGCNA_select_power(prepared, strategy = "scale_free")
modules <- WGCNA_build_modules(prepared, power, threads = 8)
sample_result <- WGCNA_module_sample(
  modules, sample_info = sample_info, time_col = "time_group"
)
membership <- WGCNA_module_membership(
  modules, MEs = sample_result$MEs_oriented
)

WGCNA_plot_soft_threshold(power, "Soft_thresholding_power.pdf")
WGCNA_export_network(
  modules, threshold = 0.15, max_edges = 5000,
  output_dir = "WGCNA_Output/Network"
)
```

提供连续或已数值编码的外部性状时，可额外运行：

```r
trait_result <- WGCNA_module_trait(
  modules, traits,
  MEs = sample_result$MEs_oriented,
  p_adjust_method = "bonferroni"
)
```

15个左右样本仅处于WGCNA建议下限附近，结果应作为探索性网络解释；客户预筛选基因得到的是候选集合内部网络，不能替代全转录组网络。对于FPKM/TPM，建议比较原值与`log2p1`结果稳定性；对于原始count，建议先在包外完成VST等方差稳定化。

## LEfSe

LEfSe 模块代码位于 `R/LEFSE.R`，从
`192.168.30.202:23099/micro_dy_gro/micro:v2.45` 的
`microbiomeMarker` 工作流提取并重构。函数不调用 Docker；在独立 R 环境安装
`microbiomeMarker`、`phyloseq`、`ggtree`、`ggplot2` 和 `dplyr` 后即可使用。

```r
colors <- c(
  Control = "#4CAF50", Model = "#3F51B5",
  `NR-HPS` = "#B30000", `R-HPS` = "#FFEB3B"
)
lefse <- LEFSE_run(otu, taxonomy, group, group_name = "Group")
markers <- LEFSE_marker_table(lefse)
LEFSE_LDA_plot(markers, colors, "LEFSE_LDA.pdf")
LEFSE_cladogram_plot(lefse, colors, "LEFSE_cladogram.pdf")
```

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
`cutoff` 是绘图前采用的筛选口径。绘图筛选严格保留 `0 < p < cutoff`，因此
`p = 0` 会被排除并给出提示，避免把数值下溢误画成无限显著；这也适用于差异
结果表按 `pvalue` 绘图时的 `p < 0.05 且 p > 0` 规则。因此可明确实现“按照
原始 P 值而不是校正后 P 值筛选和绘图”。

两种气泡图、三档 Gene count 图例和指定通路标签强调的完整示例，见
[使用 `GO_KEGG_plot()` 绘制 GO/KEGG 富集气泡图](docs/go-kegg-dotplot-tutorial.md)。

带分组色带、Gene count 气泡和基因标签的富集分析条形图示例，见
[使用 `enrichment_bar()` 绘制富集分析条形图](docs/enrichment-bar-tutorial.md)。

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

完整示例见[使用 `Heatmap_plot()` 绘制带样本分组的表达热图](docs/heatmap-tutorial.md)。

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

版本 `0.2.0` 先覆盖 GO/KEGG ORA 的计算，以及 `dotplot` 和 `barplot`。气泡图
支持 P 值、GeneRatio 或 RichFactor 横轴、三档 Count 图例，以及指定通路标签
加粗和改色。GO 三分面、圈图、结果 Excel 导出等需求将在后续版本逐步加入。

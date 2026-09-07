# GSEA 可视化与单通路批量出图

本模块借鉴返修案例 `SNCH0426260730020-0002` 的
`KEGG_GSEA_analysis.R`，并结合既有的“单通路三面板”资产整理。案例中值得复用的部分是：

- 从排序向量和 gene set 重建 running enrichment score；
- 上层 ES 曲线、中层命中条码与红–白–蓝排序色带、下层 ranked metric 紧凑展示；
- NES、P value 和 adjusted P 只标注在顶部面板；
- 可直接处理 `clusterProfiler::gseaResult`，也可处理 fgsea/普通表格场景；
- 文件名清理、批量出图和 Times New Roman 样式参数化。

没有复用案例中的客户路径、523 个筛选基因限制、人鼠 SYMBOL 大写代理或固定分组顺序。

## 函数分层

| 函数 | 适用输入 | 返回值 |
|---|---|---|
| `calculate_gsea_curve()` | 命名排序向量 + 单个 gene set | ES 曲线、hits、peak、leading edge |
| `plot_gsea()` | 已计算的 ES、hits、ranked metric | 三面板 `ggplot` |
| `plot_gsea_pathway()` | 命名排序向量 + 单个 gene set | 自动计算并绘制的 `ggplot` |
| `plot_gsea_result()` | `clusterProfiler::gseaResult` + 通路 ID/名称 | 单通路 `ggplot` |
| `plot_gsea_batch()` | `gseaResult` + 多个通路 | PDF/PNG 文件路径 |

## 1. 从 ranked vector 和 gene set 作图

这是最通用的接口，适用于 fgsea 结果、旧流程导出的向量或无法恢复 S4 对象的项目。

```r
library(BioInfoAfterSale)

set.seed(2026)
ranked_metric <- sort(rnorm(5000), decreasing = TRUE)
names(ranked_metric) <- paste0("Gene_", seq_along(ranked_metric))

# 虚构通路：偏向排名前端，并混入少量中部基因
gene_set <- unique(c(
  names(ranked_metric)[sample(1:900, 120)],
  names(ranked_metric)[sample(1800:3200, 30)]
))

style <- choose_plot_style(
  font_family = "Times New Roman",
  theme = "bw",
  dpi = 400,
  figure_width = 9.2,
  figure_height = 6.8,
  title = list(size = 18, bold = TRUE),
  axis_title = list(size = 14),
  axis_text = list(size = 12),
  data_label = list(size = 13)
)

p <- plot_gsea_pathway(
  ranked_metric = ranked_metric,
  gene_set = gene_set,
  title = "Synthetic KEGG pathway",
  statistics = c(
    NES = 2.22,
    `P value` = 6.1e-09,
    `Adjusted P` = 2.8e-06
  ),
  exponent = 1,
  colors = c("#E25659", "#335372"),
  metric_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
  style = style,
  output_file = "docs/images/GSEA.pdf",
  width = 9.2,
  height = 6.8
)

curve <- attr(p, "gsea_curve")
curve$enrichment_score
curve$peak_position
curve$leading_edge
```

结果：[GSEA.pdf](images/GSEA.pdf)

`exponent = 1` 对命中基因按 `abs(ranked_metric)` 加权；`exponent = 0` 为不加权 ES。
函数会先按分值递减排序，并拒绝重复 gene ID、无交集 gene set 以及全基因命中的无效输入。

## 2. 直接绘制 clusterProfiler GSEA 对象

```r
# gsea_result 为 clusterProfiler::GSEA() 或 gseKEGG() 的结果
p <- plot_gsea_result(
  gsea_result,
  gene_set_id = "hsa04110",
  style = style,
  output_file = "GSEA_hsa04110.pdf"
)
```

`gene_set_id` 可以是结果表中的精确 `ID` 或精确 `Description`。函数自动提取：

- `result@geneList`；
- `result@geneSets[[ID]]`；
- `NES`、`pvalue`、`p.adjust`。

批量输出客户挑选的通路：

```r
files <- plot_gsea_batch(
  result = gsea_result,
  gene_set_ids = c("hsa04110", "hsa03030", "hsa03420"),
  output_dir = "Output/GSEA_pathways",
  format = "pdf",
  style = style,
  width = 9.2,
  height = 6.8
)
```

## 3. 与 fgsea 配合

fgsea 结果表本身不保存完整 ranked vector 和 pathways，因此应把三者一起保留：

```r
fgsea_table <- as.data.frame(fgsea_result)
id <- fgsea_table$pathway[1]
row <- fgsea_table[fgsea_table$pathway == id, , drop = FALSE]

plot_gsea_pathway(
  ranked_metric = ranked_metric,
  gene_set = pathways[[id]],
  title = id,
  statistics = c(
    NES = row$NES,
    `P value` = row$pval,
    `Adjusted P` = row$padj
  ),
  style = style
)
```

新分析优先使用 `fgseaMultilevel`，不要仅为提高置换次数而给 `fgsea()` 固定 `nperm`，
否则会退回 `fgseaSimple`。

## 排序指标与解释边界

GSEA 的生物学方向由 ranked metric 的正负和排序方向决定，必须记录比较组和指标定义。

- 两组比较通常使用带方向的统计量，例如 Wald statistic、t statistic 或 signed log10 P；
- 多组 ANOVA F 值全部非负，只表示“组间差异强度”，不表示上调或下调；使用它时应采用
  正向富集语义，并避免把 NES 解释为某一组上调；
- 仅使用筛选后的差异基因会改变背景和 ES，应优先使用完整可检测基因排序；
- SYMBOL 大写不能视为可靠的人鼠同源映射，正式跨物种分析应使用明确的 ortholog 表；
- 展示的 P 值来源必须与结果表一致，不要把 raw P 和 adjusted P 混用。

## 图形参数

- `colors`：ES/正 ranked metric 与负 ranked metric 的颜色；
- `metric_colors`：中层排序色带的低值、中点和高值颜色；
- `show_metric_band = FALSE`：关闭中层色带，只保留 hits；
- `peak_position`：低层 `plot_gsea()` 中手动指定 ES 峰位置；高层接口会自动计算；
- `style`：统一控制字体、字号、边框和输出分辨率。

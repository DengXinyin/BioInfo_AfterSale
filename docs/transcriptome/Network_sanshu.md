# 转录组 PPI 与共表达网络图

`ppi_network()` 接收边表和节点注释表，使用固定随机种子的 Fruchterman–Reingold 布局。边端点必须能在节点表中完整匹配。示例网络和基因名称均为虚构数据。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "classic", dpi = 300,
  figure_width = 9, figure_height = 9,
  title = list(size = 18, bold = TRUE), legend = list(position = "right", frame = FALSE),
  panel = list(border = FALSE, major_grid = FALSE, minor_grid = FALSE)
)

set.seed(2044)
genes <- paste0("Protein", sprintf("%02d", 1:45))
edges <- data.frame(
  Node1 = sample(genes, 110, replace = TRUE),
  Node2 = sample(genes, 110, replace = TRUE),
  Score = runif(110, 0.35, 1)
)
edges <- edges[edges$Node1 != edges$Node2, ]
nodes <- data.frame(
  Gene = genes,
  Status = sample(c("Up", "Down", "NoSig"), length(genes), replace = TRUE,
                  prob = c(0.42, 0.38, 0.20))
)
ppi_network(
  edges, nodes, weight = "Score", seed = 2044, label_top = 8,
  style = style, output_file = "docs/images/transcriptome/Network_PPI_sanshu.pdf"
)
```

结果：[PPI 网络图](../images/transcriptome/Network_PPI_sanshu.pdf)

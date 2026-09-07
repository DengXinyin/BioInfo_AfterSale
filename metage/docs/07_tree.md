# 物种树与丰度

对应 `R/07_tree.R`，整理自 `bar_tree.R`。

`metage_tree_abundance()` 接收已有 `ape::phylo` 树和命名丰度向量（或特征×样本矩阵），只保留匹配的 Top N tips，返回树图和丰度条形图。`metage_tree_abundance_grob()` 将两部分排版，便于 PDF 输出。

```r
x <- metage_tree_abundance(tree, abundance, top_n = 30)
g <- metage_tree_abundance_grob(x)
metage_save_plot(g, "Output/tree_abundance.pdf", 12, 8)
```

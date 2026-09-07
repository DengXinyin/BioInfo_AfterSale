# 组成、热图与样本相关性

对应 `R/02_composition_heatmap.R`。三个函数都接收“行=物种/功能、列=样本”的丰度矩阵。

- `metage_abundance_bar()`：支持 Top N、Others、相对丰度和分组均值。
- `metage_abundance_heatmap()`：支持行 Z-score 及行/列聚类。
- `metage_sample_correlation()`：支持 Pearson/Spearman 和聚类排序。

```r
p <- metage_abundance_bar(abundance, top_n = 20, relative = TRUE)
h <- metage_abundance_heatmap(abundance, top_n = 30)
c <- metage_sample_correlation(abundance, method = "spearman")
```

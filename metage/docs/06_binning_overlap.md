# Binning 与集合关系

对应 `R/06_binning_overlap.R`。

- `metage_gc_coverage()`：输入 contig/bin 明细，列名可配置。
- `metage_bin_heatmap()`：bin×样本丰度热图。
- `metage_overlap_plot()`：2–4 个集合优先 Venn，更多集合在安装 UpSetR 时返回可绘制函数。

```r
p <- metage_gc_coverage(blob, gc = "gc", coverage = "coverage", bin = "bin")
v <- metage_overlap_plot(list(Group_A = genes_a, Group_B = genes_b))
```

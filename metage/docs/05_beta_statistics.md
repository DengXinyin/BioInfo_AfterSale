# Beta 多样性统计

对应 `R/05_beta_statistics.R`。

`metage_beta_test()` 接受 `dist` 或方阵，统一封装 ANOSIM、PERMANOVA（Adonis）与 MRPP；`permutations` 和随机种子显式记录。`metage_anosim_plot()` 返回图、模型和秩数据。

```r
d <- vegan::vegdist(t(abundance), "bray")
fit <- metage_beta_test(d, group, method = "adonis", permutations = 999)
a <- metage_anosim_plot(d, group)
```

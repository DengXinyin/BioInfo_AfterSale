# 测序质控、组装与基因长度

对应 `R/01_quality_assembly.R`，来源为 `atgc_content.R`、`error_rate.R`、`data_composition_bar.R`、`length_count.R` 和 `gene_length.R`。

- `metage_base_content_plot(data)`：首列为位置，其余数值列为 A/T/G/C/N 百分比。
- `metage_error_profile(data)`：绘制逐位点错误率。
- `metage_read_composition(counts)`：行是 reads 类别、列是样本，可绘制计数或比例。
- `metage_length_distribution(lengths)`：连续直方图；传 `breaks` 后绘制分箱堆叠图。

```r
base <- data.frame(Position = 1:100, A = runif(100, 20, 30),
                   T = runif(100, 20, 30), G = runif(100, 20, 30), C = runif(100, 20, 30))
p <- metage_base_content_plot(base)
metage_save_plot(p, "Output/base_content.pdf")
```

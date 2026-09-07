# micro_v2:v5.42 源脚本盘点

镜像最终层曾执行清空并重新复制脚本，导致流程代码仍引用的 R 文件在最终容器中不可见。
本次从该镜像的只读 OverlayFS 历史层恢复了 68 个 R 脚本（共 8319 行），再按图型合并。

## 已封装的绘图脚本

- 组成：`bar_plot.r`, `fungi_bar.R`, `func2_barplot.R`, `fun2_barplot_bugbase.R`,
  `plot_bar.r`, `heatmap.r`, `fungi_heatmap.R`, `func2_heatmap.R`, `plot_heatmap.r`,
  `tree_bar_plot.r`。
- 多样性与排序：`a_diversity.r`, `rarefaction.r`, `Rank_Abundance.r`, `Specaccum.r`,
  `PCA.r`, `PCOA.r`, `NMDS.r`, `pls_da.r`, `UPGMA.r`, `tree.r`, `dis_heatmap.r`,
  `func2_PCA.R`, `func2_PCOA.R`, `func2_NMDS.R`。
- 差异：`anova.r`, `plot_wilcoxon.r`, `func2_anova.R`, `func2_wilcoxon.R`,
  `stamp.R`, `tax_stamp.R`, `func2_stamp.R`, `fungi_stamp.R`, `plot_metastats.r`,
  `fungi_metastats.r`, `randomForest.r`, `tax_randomForest.R`, `func2_randomForest.R`,
  `fungi_randomForest.R`, `metagenomeSeq.r`, `tax_metaseq.R`, `func2_metaseq.R`,
  `fungi_metaseq.R`, `plot_ldascore.r`, `tax_LDAscore.R`, `func2_lefse.R`。
- 统计：`Anosim.r`, `Adonis.r`, `MRPP.r`, `diff_ana.r`, `tax_Anosim.R`,
  `tax_Adonis.R`, `tax_MRPP.R`, `func2_Anosim.R`, `func2_Adonis.R`, `func2_MRPP.R`,
  `fungi_Anosim.R`, `fungi_Adonis.R`, `fungi_MRPP.R`。
- 集合、网络与特殊图：`venn_upset.r`, `plot_network.r`, `zi_pi.R`, `len_Dis.r`,
  `NCM.R`, `biotype.R`, `plot-MicroPITA.R`, `plot_Cladogram.R`。

## 未直接封装为绘图函数

- `Tax4Fun2.R`, `BugBase.R`, `pre-MicroPITA.R`：分析或数据准备，不是独立绘图脚本。
- `plot_Cladogram.R`：仓库已有 `LEFSE_cladogram_plot()`，避免重复维护。
- `biotype.R`, `plot-MicroPITA.R`：高度绑定原流程中间文件；其通用排序散点图意图由
  `micro_ordination()` 覆盖。
- `fungi_Adonis.R`, `tax_Adonis.R` 等重复版本：统一由 `micro_beta_test()` 接受同一输入。

原始脚本只用于追溯和重构，没有复制进包目录，以免保留硬编码的个人路径、旧环境
`.libPaths()`、自动写文件和 HTML/PDF 副作用。

# metage:v2.87 源脚本盘点

盘点位置：`/root/microbiome/microbiome/metage_megahit`。最终镜像中共发现 39 个 R 脚本、4,482 行。原始脚本仅用于追溯和重构，没有复制进仓库，因其包含固定流程目录、自动写文件、旧字体和 HTML/PDF 副作用。

## 脚本到函数的映射

| 功能域 | 镜像原脚本 | 整理后的主要函数 |
|---|---|---|
| 测序质控 | `atgc_content.R`, `error_rate.R`, `data_composition_bar.R` | `metage_base_content_plot()`, `metage_error_profile()`, `metage_read_composition()` |
| 组装与基因长度 | `length_count.R`, `gene_length.R` | `metage_length_distribution()` |
| 组成与聚类 | `tax_bar_plot.R`, `func_barplot.R`, `tax_heatmap.R`, `func_heatmap.R`, `sample.corr_heatmap.R` | `metage_abundance_bar()`, `metage_abundance_heatmap()`, `metage_sample_correlation()` |
| 物种树+丰度 | `bar_tree.R` | `metage_tree_abundance()`, `metage_tree_abundance_grob()` |
| Alpha 多样性 | `alpha_diver.R` | `metage_alpha_indices()`, `metage_alpha_boxplot()` |
| 排序分析 | `tax_PCA.R`, `tax_PCoA.R`, `tax_NMDS.R`, `func_PCA.R`, `func_PCOA.R`, `func_NMDS.R` | `metage_ordination()`；支持优先复用已有坐标 |
| 差异箱线图 | `tax_anova.R`, `tax_wilcoxon.R`, `func_anova.R`, `func_wilcoxon.R` | `metage_diff_boxplot()` |
| STAMP | `tax_stamp.R`, `func_stamp.R` | `metage_stamp_plot()`，直接读取既有差异结果 |
| 随机森林 | `tax_randomForest.R`, `func_randomForest.R` | `metage_random_forest_importance()`，绘制既有 importance 表 |
| metagenomeSeq | `tax_metaseq.R`, `func_metaseq.R` | `metage_significant_heatmap()`；统计结果与绘图分离 |
| LEfSe | `tax_LDAscore.R`, `func_lefse.R` | `metage_lda_score_plot()` |
| Beta 统计 | `tax_Adonis.R`, `tax_Anosim.R`, `tax_MRPP.R`, `func_Adonis.R`, `func_Anosim.R`, `func_MRPP.R` | `metage_beta_test()`, `metage_anosim_plot()` |
| binning | `GC-cov.R`, `bins_heatmap.R` | `metage_gc_coverage()`, `metage_bin_heatmap()` |
| 集合关系 | `venn_flower.R` | `metage_overlap_plot()`；2–4 集合为 Venn，更多集合为 UpSet |

## 合并说明

- `tax_*` 与 `func_*` 的输入目录和文件名不同，但核心图层一致，因此按统一矩阵/结果表接口合并。
- Adonis 与 MRPP 原脚本只输出统计表，并非独立绘图；保留为模型函数，不伪装成图。
- metagenomeSeq 与 randomForest 的模型拟合不应在“重绘”时默认发生，因此函数消费既有结果；这避免改动统计口径。
- 镜像只有 LEfSe LDA 条形图 R 脚本，没有 R cladogram 绘图脚本；cladogram 可复用包根目录已有的 `LEFSE_cladogram_plot()`。

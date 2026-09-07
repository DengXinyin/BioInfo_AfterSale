# micro_v2 可视化函数库

本目录将 Docker 镜像
`192.168.30.202:23099/micro_dy_gro/micro_v2:v5.42` 中的微生物可视化脚本
整理为可直接 `source()` 的 R 函数。镜像 ID 为
`sha256:4fc6e9a49a3760d36cebd8adef2070ad680280f0e0fc27b316a1706b1a2ca606`。

## 目录

```text
micro/
├── R/          # 按功能域整理的 R 函数
├── docs/       # 与 R 文件同名的中文使用文档
├── images/     # 由虚构数据实际生成的 PNG
├── scripts/    # 一键生成示例图
└── tests/      # 最小功能测试
```

## 快速使用

```r
source("micro/R/utils.R")
source("micro/R/01_composition.R")

# 行为物种/功能，列为样本
p <- micro_abundance_bar(abundance_matrix, group = sample_group)
micro_save_plot(p, "composition.png", width = 9, height = 6)
```

全部示例图可在仓库根目录运行：

```bash
/root/anaconda3/envs/r/bin/Rscript micro/scripts/generate_examples.R
```

## 整理原则

- 输入统一为“行=特征、列=样本”的数值矩阵；分组向量建议以样本名命名。
- tax、func2、fungi 中同图型的重复脚本合并为同一个函数，不保留项目绝对路径。
- 函数返回图对象以及必要的统计表/模型，不在函数内部写死输出目录。
- 默认字体使用可移植的 `sans`；正式返修时可传 `font_family` 覆盖。
- 示例数据均为随机虚构数据，不含客户数据。

## 原脚本到函数的映射

| 功能域 | 镜像原脚本 | 整理后的主要函数 |
|---|---|---|
| 群落组成 | `bar_plot.r`, `heatmap.r`, `tree_bar_plot.r`, `fungi_bar.R`, `fungi_heatmap.R`, `func2_barplot.R`, `func2_heatmap.R`, `plot_bar.r`, `plot_heatmap.r` | `micro_abundance_bar()`, `micro_abundance_heatmap()`, `micro_clustered_abundance_bar()` |
| Alpha/抽样 | `a_diversity.r`, `rarefaction.r`, `Rank_Abundance.r`, `Specaccum.r` | `micro_alpha_boxplot()`, `micro_rarefaction_curve()`, `micro_rank_abundance()`, `micro_species_accumulation()` |
| 排序/距离 | `PCA.r`, `PCOA.r`, `NMDS.r`, `pls_da.r`, `UPGMA.r`, `tree.r`, `dis_heatmap.r` 及 func2 版本 | `micro_ordination()`, `micro_distance_heatmap()` |
| 差异分析 | `anova.r`, `plot_wilcoxon.r`, `stamp.R`, `tax_stamp.R`, `func2_stamp.R`, `fungi_stamp.R`, `randomForest.r` 及 tax/func2/fungi 版本、`plot_ldascore.r`, `tax_LDAscore.R` | `micro_diff_boxplot()`, `micro_stamp_plot()`, `micro_random_forest_plot()`, `micro_lda_score_plot()` |
| Beta 统计 | `Anosim.r`, `Adonis.r`, `MRPP.r`, `diff_ana.r` 及 tax/func2/fungi 版本 | `micro_beta_test()`, `micro_anosim_plot()` |
| 集合/网络 | `venn_upset.r`, `plot_network.r`, `zi_pi.R` | `micro_overlap_plot()`, `micro_network_plot()`, `micro_zipi_plot()` |
| 特殊图 | `len_Dis.r`, `NCM.R` | `micro_length_distribution()`, `micro_ncm_fit()`, `micro_ncm_plot()` |

`Tax4Fun2.R`、`BugBase.R`、`pre-MicroPITA.R` 等只负责分析或数据准备的脚本未冒充绘图函数；
LEfSe cladogram 已由包根目录现有的 `LEFSE_cladogram_plot()` 覆盖。完整盘点见
[SOURCE_INVENTORY.md](SOURCE_INVENTORY.md)。

# 转录组可视化索引

本目录按功能组织转录组教程。归档和指定流程镜像中的数据图已按下表映射；通用图形复用仓库现有函数，新能力由文件名以 `_sanshu.R` 结尾的源码提供。

| 报告内容 | 教程 | 主要函数 |
|---|---|---|
| 测序错误率、碱基含量、基因组区域、覆盖度、饱和度 | [QC](QC_sanshu.md) | `read_qc_error()`、`read_qc_base_content()`、`mapping_region_pie()`、`gene_body_coverage()`、`junction_saturation()` |
| 注释花瓣、KOG/GO/KEGG/TF 柱图、NR 饼图、新转录本注释 | [注释汇总](Annotation_sanshu.md) | `annotation_flower()`、`annotation_bar()`、`annotation_pie()` |
| FPKM 箱线/密度、相关性、PCA | [表达量](Expression_sanshu.md) | `plot_violin_box()`、`plot_distribution()`、`sample_correlation_heatmap()`、`plot_pca()` |
| DEG 数量、Venn/花瓣、UpSet、MA、火山、热图、趋势 | [差异表达](DEG_sanshu.md) | `deg_count_bar()`、`deg_flower()`、`deg_upset()`、`deg_ma()`、`plot_volcano()`、`Heatmap_plot()`、`expression_trend()` |
| GO/KEGG 蝴蝶、条形、气泡、DAG、圈图 | [富集分析](Enrichment_sanshu.md) | `enrichment_butterfly()`、`GO_KEGG_plot()`、`go_dag()`、`enrichment_circle()` |
| GSEA 运行分数、GSVA 热图与差异通路 | [GSEA](../GSEA.md)、[通用热图](../Heatmap.md)、[注释柱图](Annotation_sanshu.md) | `plot_gsea()`、`Heatmap_plot()`、`annotation_bar()` |
| PPI 和共表达网络 | [网络](Network_sanshu.md) | `ppi_network()` |
| 可变剪接、DEU、SNP/INDEL | [基因结构与变异](GeneStructure_Variant_sanshu.md) | `annotation_bar()`、`deu_exon_expression()`、`annotation_pie()` |
| 样本树、软阈值、网络、模块数量、模块性状、热图、MM–GS | [WGCNA](WGCNAPlots_sanshu.md) | `wgcna_sample_tree()`、`WGCNA_plot_soft_threshold()`、`ppi_network()`、`wgcna_module_sizes()`、`wgcna_module_trait_heatmap()`、`Heatmap_plot()`、`wgcna_mm_gs()` |

概念流程图、FASTQ 格式示意图和可变剪接类型示意图不是数据可视化函数。KEGG 原生 pathway map 缺少版本匹配的 KGML/底图时不可重建；大型 WGCNA TOM 网络缺少边表时也不可从报告图片反推。

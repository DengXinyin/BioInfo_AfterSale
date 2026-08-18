# WGCNA实际数据示例教程

本教程使用此前项目中的935个筛选基因和15个时间序列样本，演示无外部表型WGCNA。完整可执行版本位于`inst/examples/WGCNA_example_tutorial.R`。

## 1. 数据与环境

示例输入：

```text
/home/xydeng/SNKH042726062901_余杰辉_转录生信分析售后/
├── Input/selected_gene_expression_fpkm.tsv  # 935基因×15样本
└── Input/sampleinfo.tsv                     # 样本、时间点、小时、重复
```

R包依赖：

```r
# 非base直接依赖；只需在环境构建时安装一次
install.packages("WGCNA")

# 已安装完整R包时
library(BioInfoAfterSale)

# 在v1.6.7镜像中只验证WGCNA模块时，可直接加载源码
source("/home/xydeng/Useful_Docker_Images/BioInfo_AfterSale/R/WGCNA.R")
```

函数内部已经使用`requireNamespace("WGCNA")`检查依赖，并通过`WGCNA::函数`调用。不会在正式分析过程中联网安装R包。

## 2. 读取数据

```r
project_dir <- "/home/xydeng/SNKH042726062901_余杰辉_转录生信分析售后"
expression <- read.delim(
  file.path(project_dir, "Input/selected_gene_expression_fpkm.tsv"),
  row.names = 1,
  check.names = FALSE
)
sample_info <- read.delim(
  file.path(project_dir, "Input/sampleinfo.tsv"),
  check.names = FALSE
)

stopifnot(nrow(expression) == 935, ncol(expression) == 15)
stopifnot(setequal(colnames(expression), sample_info$Sample_ID))
```

## 3. 一键复刻v1.6.7镜像逻辑

客户要求附件中全部基因纳入，因此关闭零值比例过滤，并令`goodSamplesGenes()`不允许静默删除数据。`strategy = "image_rule"`按镜像规则在15个样本时选择unsigned power 9。

```r
output_dir <- "/tmp/BioInfoAfterSale_WGCNA_Tutorial"

result <- WGCNA_run(
  expression,
  sample_info = sample_info,
  time_col = "time_group",
  prepare_args = list(
    transform = "none",
    keep_all_genes = TRUE,
    remove_bad = FALSE
  ),
  power_args = list(
    network_type = "unsigned",
    strategy = "image_rule",
    fit_cutoff = 0.9
  ),
  module_args = list(
    deep_split = 2,
    min_module_size = 30,
    merge_cut_height = 0.25,
    reassign_threshold = 0,
    pam_respects_dendro = FALSE,
    max_block_size = 1000,
    threads = 8,
    seed = 123
  ),
  output_dir = output_dir
)
```

本数据应得到power 9、7个模块（含grey）和935个模块分配记录。镜像复刻结果的模块大小为turquoise 316、blue 238、brown 132、yellow 87、green 73、red 72、grey 17。

## 4. 查看模块—样本关系

没有外部表型时，不构造“一个样本为1、其他样本为0”的伪性状。使用模块eigengene（ME）表达模块在各样本中的整体变化：

```r
ME <- result$module_sample$MEs_oriented
ME_zscore <- result$module_sample$MEs_zscore
module_sample_long <- result$module_sample$sample_values
module_time_summary <- result$module_sample$time_summary

head(module_sample_long)
head(module_time_summary)
```

`MEs_raw`保留WGCNA原始主成分方向；`MEs_oriented`被统一定向为与模块平均表达正相关。实际翻转情况见：

```r
result$module_sample$orientation
```

如需复用包内热图函数：

```r
Heatmap_plot(
  t(ME_zscore),
  group = setNames(sample_info$time_group, sample_info$Sample_ID),
  scale = "none",
  output_file = file.path(output_dir, "module_sample_ME_heatmap.pdf"),
  figure_width = 12,
  figure_height = 7
)
```

## 5. soft power诊断

```r
result$power$power
result$power$reason
result$power$fit_indices

WGCNA_plot_soft_threshold(
  result$power,
  file.path(output_dir, "Soft_thresholding_power.pdf")
)
```

`image_rule`用于复刻云流程。若新项目希望按诊断选择power，改用：

```r
power_diagnostic <- WGCNA_select_power(
  result$prepared,
  network_type = "unsigned",
  strategy = "scale_free",
  fit_cutoff = 0.9
)
```

## 6. 模块基因和hub gene

```r
module_gene <- result$modules$module_gene
assigned_kME <- result$membership$assigned

# blue模块按|kME|排序的前20个候选hub gene
blue_hub <- subset(assigned_kME, module == "blue")
blue_hub <- blue_hub[order(blue_hub$rank_in_module), ]
head(blue_hub, 20)
```

`grey`表示未稳定归入共表达模块，不建议作为功能模块解释。

## 7. 导出TOM网络

以下示例只导出blue模块最强的5000条边，避免网络图过密：

```r
network <- WGCNA_export_network(
  result$modules,
  module = "blue",
  threshold = 0.15,
  max_module_genes = 2000,
  max_edges = 5000,
  output_dir = file.path(output_dir, "Network")
)

head(network$blue$nodes)
head(network$blue$edges)
```

该网络是表达数据计算出的TOM共表达网络，不是STRING或其他PPI网络。

## 8. 有真实外部性状时

仅在性状具有真实实验含义时计算模块—性状相关。例如连续性状表必须以样本为行，并与表达矩阵样本名匹配：

```r
# traits <- data.frame(
#   Sample_ID = sample_info$Sample_ID,
#   disease_score = measured_disease_score
# )
# trait_result <- WGCNA_module_trait(
#   result$modules,
#   traits,
#   MEs = result$module_sample$MEs_oriented,
#   p_adjust_method = "bonferroni"
# )
# trait_result$table
```

## 9. FPKM尺度敏感性分析

当前一键示例使用原始FPKM以复刻镜像。FPKM跨度较大时，建议另跑`log2(FPKM+1)`版本，并比较模块数、grey比例及模块成员稳定性：

```r
log_result <- WGCNA_run(
  expression,
  sample_info = sample_info,
  time_col = "time_group",
  prepare_args = list(
    transform = "log2p1",
    keep_all_genes = TRUE,
    remove_bad = FALSE
  ),
  power_args = list(
    network_type = "unsigned",
    strategy = "scale_free",
    fit_cutoff = 0.9
  ),
  module_args = list(threads = 8),
  output_dir = file.path(output_dir, "log2p1_sensitivity")
)
```

不能因为镜像有固定默认值就忽略数据尺度和soft power诊断。最终报告应明确表达值变换、network类型、power及其选择依据。

## 10. 自动生成的核心文件

`WGCNA_write_results()`或`WGCNA_run(output_dir=...)`会生成：

- `expression_qc.csv`
- `soft_threshold_fit.csv`
- `modulegene.csv`
- `module_summary.csv`
- `module_eigengene_raw.csv`
- `module_eigengene_oriented.csv`
- `module_sample_eigengene_long.csv`
- `module_time_summary.csv`
- `gene_module_correlation.csv`
- `assigned_module_kME.csv`
- `RUN_MANIFEST.csv`

建议同时保存调用脚本、标准错误日志和`sessionInfo()`，确保售后项目可追溯。

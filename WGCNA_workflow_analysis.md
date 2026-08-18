# WGCNA流程分析与R包接口设计

## 1. 参考环境

- 镜像：`192.168.30.202:23099/wgcna/wgcna:v1.6.7`
- 镜像脚本：`/script/wgcna2.r`
- R：4.1.2
- WGCNA：1.73
- `v1.6.7`与`v1.6.4`中的`/script/wgcna2.r` SHA256均为
  `ca854c02132d7a8daa16b1c591c631bb79d182d3a9718cee36b201d3653c509d`。

新模块不在函数内部启动Docker，而是在安装了WGCNA的R环境中直接运行。镜像用于确定算法、默认参数和兼容性基线。

## 2. 镜像原流程

镜像脚本依次完成：

1. 读取基因×样本表达矩阵、样本性状表和样本分组表；
2. 删除零值比例超过30%的基因；
3. 转置为样本×基因并运行`goodSamplesGenes()`；
4. 样本聚类；
5. 在1–10及12、14、16、18、20中评估soft power；
6. 不使用诊断结果自动选power，而是按样本数固定power；
7. 用`blockwiseModules()`构建模块；
8. 输出模块基因、树图、模块大小和抽样TOM热图；
9. 计算模块—性状、gene significance、module membership；
10. 按模块重算TOM并导出Cytoscape网络；
11. 输出模块基因列表，供GO/KEGG富集使用。

镜像关键默认值如下：

| 参数 | 默认值 |
|---|---|
| network/TOM | unsigned |
| deepSplit | 2 |
| minModuleSize | 30 |
| mergeCutHeight | 0.25 |
| reassignThreshold | 0 |
| pamRespectsDendro | FALSE |
| maxBlockSize | 1000 |
| 随机种子 | 123 |
| TOM导出阈值 | 0.15 |

样本数soft power规则：

| 样本数 | unsigned | signed/signed hybrid |
|---|---:|---:|
| <20 | 9 | 18 |
| 20–29 | 8 | 16 |
| 30–39 | 7 | 14 |
| ≥40 | 6 | 12 |

## 3. R包拆分

| 流程步骤 | 函数 | 主要返回值 |
|---|---|---|
| 输入/QC | `WGCNA_prepare_expression()` | `datExpr`、QC表、剔除ID、样本树 |
| soft power | `WGCNA_select_power()` | fit indices、power、选择原因 |
| 模块构建 | `WGCNA_build_modules()` | 原始net、modulegene、ME、模块汇总 |
| 无表型样本关系 | `WGCNA_module_sample()` | 原始/定向ME、长表、分组均值±SD |
| 模块—性状 | `WGCNA_module_trait()` | r、P、校正P及长表 |
| gene—module | `WGCNA_module_membership()` | 全部kME及所属模块排名 |
| TOM网络 | `WGCNA_export_network()` | nodes/edges及可选文件 |
| 写结果 | `WGCNA_write_results()` | 标准CSV结果集 |
| soft power图 | `WGCNA_plot_soft_threshold()` | PDF/PNG诊断图 |
| 一键入口 | `WGCNA_run()` | 全部中间对象和结果 |

所有新增公开函数均以`WGCNA_`开头。计算函数返回对象而不是只写文件，便于后续R包继续绘图、筛选hub gene或衔接富集分析。

## 4. 相对镜像的通用化改造

### 4.1 soft power策略显式化

- `strategy = "image_rule"`：严格复刻镜像的样本数规则；
- `strategy = "scale_free"`：选择首个达到R²阈值且斜率为负的候选；没有候选时回退镜像规则；
- `strategy = "manual"`：由调用者指定，并在返回对象中保留选择原因。

不能只保存soft power图而不记录实际power及选取依据。

### 4.2 表达值变换显式化

- `transform = "none"`：复刻镜像；
- `transform = "log2p1"`：适合需要压缩FPKM/TPM长尾时；
- 原始count建议先在包外进行VST等处理，再传给WGCNA。

表达变换会改变相关结构。正式项目应根据数据来源选择，并在必要时比较模块稳定性。

### 4.3 所有客户基因纳入

镜像默认删除零值比例超过30%的基因。客户明确要求所有基因纳入时，应使用：

```r
WGCNA_prepare_expression(
  expression,
  keep_all_genes = TRUE,
  remove_bad = FALSE
)
```

`keep_all_genes`只关闭零比例规则；`goodSamplesGenes()`仍执行。若数据无法通过且`remove_bad = FALSE`，函数停止而不是静默违背客户要求。

### 4.4 无外部表型

不把单个样本编码为一个1、其余样本为0后计算所谓模块—样本相关P值。这类结果本质上是标准化ME值的变形，缺乏可解释的组间统计基础。

无表型时使用：

- 每个模块在每个样本中的ME；
- 模块×样本ME热图；
- 有重复分组时的描述性ME均值±SD；
- 模块内kME及hub gene排序。

若有真正的连续性状或已明确编码的实验组变量，再使用`WGCNA_module_trait()`。

### 4.5 ME符号

主成分符号可整体乘以-1。`WGCNA_module_sample()`同时保留原始ME，并可将展示ME定向为与模块平均表达正相关；翻转倍数和相关系数均保存在`orientation`表中。

## 5. 标准调用

```r
result <- WGCNA_run(
  expression,
  sample_info = sample_info,
  time_col = "time_group",
  prepare_args = list(
    transform = "log2p1",
    keep_all_genes = TRUE,
    remove_bad = FALSE
  ),
  power_args = list(
    strategy = "scale_free",
    network_type = "unsigned",
    fit_cutoff = 0.9
  ),
  module_args = list(
    deep_split = 2,
    min_module_size = 30,
    merge_cut_height = 0.25,
    threads = 8
  ),
  output_dir = "WGCNA_Output"
)
```

带真实外部性状时传入`traits`；需要Cytoscape网络时设置`export_networks = TRUE`并在`network_args`中配置阈值和边数。

## 6. 解释边界

- 建议至少20个样本；约15个样本仅处于可运行下限，模块稳定性需保守解释；
- 客户预筛选基因建立的是候选集合内部网络，不能替代全转录组网络；
- `unsigned`会把强负相关也视为强连接，时间序列项目可比较`signed`；
- grey表示未稳定归入模块，不建议作为一个功能模块解释；
- TOM网络是表达共表达网络，不是STRING/PPI；
- 网络图截取top-N边只影响展示，正式边表是否截取由`max_edges`明确控制。

## 7. 已完成验证

- 在v1.6.7镜像中用935基因×15样本实际数据运行`image_rule`；
- 得到power 9、7个模块（含grey）以及316/238/132/87/73/72/17的模块大小，与原镜像分析结果一致；
- 15×7 ME、935行assigned kME均通过维度和缺失值检查；
- 使用合成双模块数据验证模块构建、样本汇总、模块—性状、kME、网络导出和soft power PDF。

镜像不包含本R包其余模块所需的`clusterProfiler`、`ComplexHeatmap`、`microbiomeMarker`等依赖，因此目前采用“单独source WGCNA模块并在镜像验证”的方式，而不是在该镜像内安装整个R包。

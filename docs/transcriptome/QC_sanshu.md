# 转录组测序与比对质控图

本教程覆盖双端测序错误率、碱基含量、基因组区域分布、基因体覆盖度和剪接位点饱和度。示例数据均为固定随机种子生成的虚构数据：每行分别表示碱基位置、区域计数或某一曲线观测点，不含真实项目样本及数值。示例只输出 PDF。

```r
library(BioInfoAfterSale)

dir.create("docs/images/transcriptome", recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  font_family = "Times New Roman", theme = "bw", dpi = 300,
  figure_width = 10, figure_height = 6,
  title = list(size = 18, bold = TRUE),
  axis_title = list(size = 13), axis_text = list(size = 10),
  legend_title = list(size = 12), legend_text = list(size = 10),
  legend = list(position = "right", frame = FALSE)
)

set.seed(2041)
qc <- expand.grid(
  SampleID = "RNA-A01", Stage = c("RawData", "CleanData"),
  Read = c("read1", "read2"), Position = 1:100,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
qc$ErrorRate <- pmax(0.02, 0.18 + qc$Position / 900 +
  ifelse(qc$Stage == "RawData", 0.08, 0) + rnorm(nrow(qc), 0, 0.012))
read_qc_error(
  qc, read = "Read", colors = c(RawData = "#9ECFC0", CleanData = "#65A479"),
  style = style, output_file = "docs/images/transcriptome/QC_error_sanshu.pdf"
)

base_content <- expand.grid(
  SampleID = "RNA-A01", Stage = c("RawData", "CleanData"),
  Read = c("read1", "read2"), Position = 1:100,
  Base = c("A", "T", "C", "G", "N", "GC"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
base_mean <- c(A = 25, T = 25, C = 24, G = 24, N = 0.3, GC = 48)
base_content$Ratio <- base_mean[base_content$Base] +
  1.2 * sin(base_content$Position / 13) + rnorm(nrow(base_content), 0, 0.25)
read_qc_base_content(
  base_content, read = "Read", style = style,
  output_file = "docs/images/transcriptome/QC_base_content_sanshu.pdf"
)

region <- data.frame(
  Sample = "RNA-A01", Region = c("exonic", "intronic", "intergenic"),
  Percentage = c(71, 17, 12)
)
mapping_region_pie(
  region, sample = "Sample",
  colors = c(exonic = "#8DD3C7", intronic = "#FFFFB3", intergenic = "#BEBADA"),
  style = style, output_file = "docs/images/transcriptome/QC_mapping_region_sanshu.pdf"
)

coverage <- data.frame(
  SampleID = rep(c("RNA-A01", "RNA-B01"), each = 101),
  GeneBodyPercent = rep(0:100, 2)
)
coverage$Coverage <- c(
  0.42 + 0.5 * (coverage$GeneBodyPercent[1:101] / 100)^0.7,
  0.48 + 0.42 * (coverage$GeneBodyPercent[102:202] / 100)^0.8
)
gene_body_coverage(
  coverage, style = style,
  output_file = "docs/images/transcriptome/QC_gene_body_coverage_sanshu.pdf"
)

saturation <- expand.grid(
  SampleID = "RNA-A01", PercentReads = seq(5, 100, 5),
  JunctionType = c("All junctions", "known junctions", "novel junctions"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
limit <- c("All junctions" = 82, "known junctions" = 61, "novel junctions" = 22)
saturation$JunctionCount <- limit[saturation$JunctionType] *
  (1 - exp(-saturation$PercentReads / 25))
junction_saturation(
  saturation, style = style,
  output_file = "docs/images/transcriptome/QC_junction_saturation_sanshu.pdf"
)
```

结果：

- [错误率](../images/transcriptome/QC_error_sanshu.pdf)
- [碱基含量](../images/transcriptome/QC_base_content_sanshu.pdf)
- [基因组区域分布](../images/transcriptome/QC_mapping_region_sanshu.pdf)
- [基因体覆盖度](../images/transcriptome/QC_gene_body_coverage_sanshu.pdf)
- [剪接位点饱和度](../images/transcriptome/QC_junction_saturation_sanshu.pdf)

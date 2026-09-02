# Polar volcano preview using the real PBMC3K differential-expression result.
input_tsv <- "docs/images/pbmc3k-cd14mono-vs-naivecd4-gene-umap.tsv"
output_png <- "docs/images/pbmc3k-cd14mono-vs-naivecd4-polar-volcano.png"
output_pdf <- "docs/images/pbmc3k-cd14mono-vs-naivecd4-polar-volcano.pdf"

data <- utils::read.delim(input_tsv, check.names = FALSE)
source("R/utils.R")
source("R/choose_plot_style.R")
source("R/volcano_plot.R")
source("R/volcano_polar_plot.R")
stopifnot(capabilities("cairo"))

plot <- volcano_polar_plot(
  result = data,
  gene_column = "gene",
  log2fc_column = "avg_log2FC",
  padj_column = "p_val_adj",
  log2fc_cutoff = 1,
  padj_cutoff = 0.05,
  title = "Polar volcano: CD14+ Mono vs Naive CD4 T",
  output_file = output_png,
  figure_width = 8,
  figure_height = 8,
  dpi = 300
)
ggplot2::ggsave(
  output_pdf, plot = plot, width = 8, height = 8,
  units = "in", device = grDevices::cairo_pdf
)

counts <- table(plot$data$.Group)
stopifnot(
  all(is.finite(plot$data$.polar_x)),
  all(is.finite(plot$data$.polar_y)),
  all(plot$data$.padj > 0)
)
cat("genes_plotted=", nrow(plot$data), "\n", sep = "")
cat("Up=", unname(counts[["Up"]]), "\n", sep = "")
cat("Down=", unname(counts[["Down"]]), "\n", sep = "")
cat("Not_significant=", unname(counts[["Not significant"]]), "\n", sep = "")
cat("fc_limit=", signif(attr(plot, "fc_limit"), 4), "\n", sep = "")
cat("png=", normalizePath(output_png), "\n", sep = "")
cat("pdf=", normalizePath(output_pdf), "\n", sep = "")
cat("polar_volcano_check=PASS\n")

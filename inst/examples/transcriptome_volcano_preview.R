# Preview the package volcano plot with reproducible transcriptome-like data.
set.seed(20260825)

n_genes <- 12000L
result <- data.frame(
  gene = sprintf("Gene_%05d", seq_len(n_genes)),
  log2FoldChange = rnorm(n_genes, mean = 0, sd = 0.72)
)

# Add asymmetric biological-effect tails to resemble a typical RNA-seq result.
up_ids <- sample(seq_len(n_genes), 620L)
down_ids <- sample(setdiff(seq_len(n_genes), up_ids), 470L)
result$log2FoldChange[up_ids] <- rnorm(length(up_ids), 1.85, 0.55)
result$log2FoldChange[down_ids] <- rnorm(length(down_ids), -1.75, 0.50)

signal <- pmax(abs(result$log2FoldChange) - 0.35, 0)
z_score <- signal * runif(n_genes, 1.8, 3.8) + rexp(n_genes, rate = 2.2)
raw_p <- pmin(runif(n_genes), 2 * stats::pnorm(-z_score))
effect_ids <- c(up_ids, down_ids)
raw_p[effect_ids] <- pmin(
  raw_p[effect_ids],
  10^runif(length(effect_ids), -11, -2.2)
)
result$padj <- p.adjust(raw_p, method = "BH")

# Mimic numerical underflow sometimes present in differential-expression files.
# volcano_plot() must remove these rows before applying -log10(padj).
zero_padj_ids <- sample(seq_len(n_genes), 36L)
result$padj[zero_padj_ids] <- 0

source("R/utils.R")
source("R/choose_plot_style.R")
source("R/volcano_normal.R")
source("R/volcano_candy_plot.R")
source("R/volcano_polar_plot.R")

stopifnot(capabilities("cairo"))

preview_png <- "docs/images/transcriptome-volcano-preview.png"
preview_pdf <- "docs/images/transcriptome-volcano-preview.pdf"
candy_png <- "docs/images/transcriptome-volcano-candy-preview.png"
candy_pdf <- "docs/images/transcriptome-volcano-candy-preview.pdf"
polar_png <- "docs/images/transcriptome-volcano-polar-preview.png"
polar_pdf <- "docs/images/transcriptome-volcano-polar-preview.pdf"

plot <- volcano_plot(
  result = result,
  log2fc_cutoff = 1,
  padj_cutoff = 0.05,
  title = "Transcriptome differential expression",
  font_family = "sans",
  point_size = 1.65,
  point_alpha = 0.68,
  output_file = preview_png,
  figure_width = 8,
  figure_height = 7,
  dpi = 300
)

ggplot2::ggsave(
  preview_pdf, plot = plot, width = 8, height = 7,
  units = "in", device = grDevices::cairo_pdf
)

candy_plot <- volcano_candy_plot(
  result = result,
  gene_column = "gene",
  log2fc_cutoff = 1,
  padj_cutoff = 0.05,
  n_neighbors = 40,
  min_dist = 0.16,
  target_weight = 0.68,
  candy_strength = 0.88,
  label_top = 0,
  point_size = 0.9,
  point_alpha = 0.58,
  output_file = candy_png,
  figure_width = 9,
  figure_height = 6.5,
  dpi = 300
)
ggplot2::ggsave(
  candy_pdf, plot = candy_plot, width = 9, height = 6.5,
  units = "in", device = grDevices::cairo_pdf
)

polar_plot <- volcano_polar_plot(
  result = result,
  gene_column = "gene",
  log2fc_cutoff = 1,
  padj_cutoff = 0.05,
  point_size = 0.9,
  point_alpha = 0.58,
  title = "Polar volcano: transcriptome differential expression",
  output_file = polar_png,
  figure_width = 8,
  figure_height = 8,
  dpi = 300
)
ggplot2::ggsave(
  polar_pdf, plot = polar_plot, width = 8, height = 8,
  units = "in", device = grDevices::cairo_pdf
)

counts <- table(plot$data$.Group)
cat("input_rows=", nrow(result), "\n", sep = "")
cat("valid_rows=", nrow(plot$data), "\n", sep = "")
cat("Up=", unname(counts[["Up"]]), "\n", sep = "")
cat("Down=", unname(counts[["Down"]]), "\n", sep = "")
cat("Not_significant=", unname(counts[["Not significant"]]), "\n", sep = "")
cat("png=", normalizePath(preview_png), "\n", sep = "")
cat("pdf=", normalizePath(preview_pdf), "\n", sep = "")
cat("candy_png=", normalizePath(candy_png), "\n", sep = "")
cat("candy_pdf=", normalizePath(candy_pdf), "\n", sep = "")
cat("polar_png=", normalizePath(polar_png), "\n", sep = "")
cat("polar_pdf=", normalizePath(polar_pdf), "\n", sep = "")

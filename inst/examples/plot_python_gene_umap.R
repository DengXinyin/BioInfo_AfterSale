# Plot Python-computed gene UMAP coordinates in R.
input_tsv <- "docs/images/pbmc3k-python-gene-umap.tsv"
output_png <- "docs/images/pbmc3k-python-r-candy-volcano.png"
output_pdf <- "docs/images/pbmc3k-python-r-candy-volcano.pdf"

data <- utils::read.delim(input_tsv, check.names = FALSE)
source("R/utils.R")
source("R/choose_plot_style.R")
source("R/volcano_plot.R")
source("R/volcano_candy_plot.R")
stopifnot(capabilities("cairo"))

plot <- volcano_candy_plot(
  result = data,
  gene_column = "gene",
  log2fc_column = "avg_log2FC",
  padj_column = "p_val_adj",
  log2fc_cutoff = 1,
  padj_cutoff = 0.05,
  umap_columns = c("UMAP_1", "UMAP_2"),
  label_top = 0,
  point_size = 0.9,
  point_alpha = 0.58,
  title = "CD14+ Mono vs Naive CD4 T",
  output_file = output_png,
  figure_width = 9,
  figure_height = 6.5,
  dpi = 300
)
ggplot2::ggsave(
  output_pdf, plot = plot, width = 9, height = 6.5,
  units = "in", device = grDevices::cairo_pdf
)

counts <- table(plot$data$.Group)
centres <- stats::aggregate(
  cbind(.UMAP_1, .UMAP_2) ~ .Group, data = plot$data, FUN = mean
)
stopifnot(
  all(is.finite(plot$data$.UMAP_1)),
  all(is.finite(plot$data$.UMAP_2)),
  centres$.UMAP_1[centres$.Group == "Down"] <
    centres$.UMAP_1[centres$.Group == "Not significant"],
  centres$.UMAP_1[centres$.Group == "Not significant"] <
    centres$.UMAP_1[centres$.Group == "Up"]
)
cat("Up=", unname(counts[["Up"]]), "\n", sep = "")
cat("Down=", unname(counts[["Down"]]), "\n", sep = "")
cat("Not_significant=", unname(counts[["Not significant"]]), "\n", sep = "")
cat("png=", normalizePath(output_png), "\n", sep = "")
cat("pdf=", normalizePath(output_pdf), "\n", sep = "")
cat("python_to_r_check=PASS\n")

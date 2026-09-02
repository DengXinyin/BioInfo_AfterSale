# Real-data preview: CD14+ monocytes versus naive CD4 T cells in PBMC3K.
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Package 'Seurat' is required for this example.")
}

suppressPackageStartupMessages(library(Seurat))
set.seed(20260825)

input_rds <- "/Users/hsinyinteng/Multi_Omics/scRNA_seq/Seurat/pbmc3k_final.rds"
object <- readRDS(input_rds)

ident_1 <- "CD14+ Mono"
ident_2 <- "Naive CD4 T"
cells_per_group <- 400L
cells_1 <- sample(WhichCells(object, idents = ident_1), cells_per_group)
cells_2 <- sample(WhichCells(object, idents = ident_2), cells_per_group)
subset_object <- subset(object, cells = c(cells_1, cells_2))

de <- FindMarkers(
  subset_object,
  ident.1 = ident_1,
  ident.2 = ident_2,
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0,
  only.pos = FALSE,
  verbose = FALSE
)
de$gene <- rownames(de)

source("R/utils.R")
source("R/choose_plot_style.R")
source("R/volcano_plot.R")
source("R/volcano_candy_plot.R")
stopifnot(capabilities("cairo"))

output_png <- "docs/images/pbmc3k-cd14mono-vs-naivecd4-candy-volcano.png"
output_pdf <- "docs/images/pbmc3k-cd14mono-vs-naivecd4-candy-volcano.pdf"
output_tsv <- "docs/images/pbmc3k-cd14mono-vs-naivecd4-gene-umap.tsv"

plot <- volcano_candy_plot(
  result = de,
  gene_column = "gene",
  log2fc_column = "avg_log2FC",
  padj_column = "p_val_adj",
  log2fc_cutoff = 0.25,
  padj_cutoff = 0.05,
  n_neighbors = 40,
  min_dist = 0.16,
  target_weight = 0.68,
  candy_strength = 0.88,
  seed = 20260825,
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

coordinate_table <- plot$data[, c(
  "gene", "avg_log2FC", "p_val_adj", ".minus_log10_padj", ".Group",
  ".UMAP_raw_1", ".UMAP_raw_2", ".UMAP_1", ".UMAP_2"
)]
names(coordinate_table) <- c(
  "gene", "avg_log2FC", "p_val_adj", "minus_log10_padj", "group",
  "UMAP_raw_1", "UMAP_raw_2", "UMAP_1", "UMAP_2"
)
utils::write.table(
  coordinate_table, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE
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

cat("input=", normalizePath(input_rds), "\n", sep = "")
cat("comparison=", ident_1, " vs ", ident_2, "\n", sep = "")
cat("cells_per_group=", cells_per_group, "\n", sep = "")
cat("genes_tested=", nrow(de), "\n", sep = "")
cat("padj_zero=", sum(de$p_val_adj == 0, na.rm = TRUE), "\n", sep = "")
cat("Up=", unname(counts[["Up"]]), "\n", sep = "")
cat("Down=", unname(counts[["Down"]]), "\n", sep = "")
cat("Not_significant=", unname(counts[["Not significant"]]), "\n", sep = "")
cat("png=", normalizePath(output_png), "\n", sep = "")
cat("pdf=", normalizePath(output_pdf), "\n", sep = "")
cat("coordinates=", normalizePath(output_tsv), "\n", sep = "")
cat("gene_umap_check=PASS\n")

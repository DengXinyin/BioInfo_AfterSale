# Static vector-PDF view of the unsupervised three-dimensional gene UMAP.
if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
  stop("Package 'scatterplot3d' is required.")
}

input_tsv <- "docs/images/pbmc3k-expression-gene-umap-3d.tsv"
output_pdf <- "docs/images/pbmc3k-expression-gene-umap-3d.pdf"
data <- utils::read.delim(input_tsv, check.names = FALSE)

groups <- c("Not significant", "Down", "Up")
data$group <- factor(data$group, levels = groups)
data <- data[order(data$group), , drop = FALSE]
colors <- c(
  Up = grDevices::adjustcolor("#F04438", alpha.f = 0.72),
  Down = grDevices::adjustcolor("#3977B8", alpha.f = 0.78),
  `Not significant` = grDevices::adjustcolor("#B7B7B7", alpha.f = 0.34)
)
point_colors <- unname(colors[as.character(data$group)])

stopifnot(
  capabilities("cairo"),
  all(is.finite(data$UMAP_1)),
  all(is.finite(data$UMAP_2)),
  all(is.finite(data$UMAP_3))
)
grDevices::cairo_pdf(output_pdf, width = 10, height = 8, family = "sans")
graphics::par(mar = c(3.2, 3.2, 4.5, 2), bg = "white")
scatterplot3d::scatterplot3d(
  x = data$UMAP_1,
  y = data$UMAP_2,
  z = data$UMAP_3,
  color = point_colors,
  pch = 16,
  cex.symbols = 0.42,
  angle = 48,
  scale.y = 0.82,
  grid = TRUE,
  box = FALSE,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  zlab = "UMAP 3",
  main = "3D expression-based gene UMAP",
  sub = "CD14+ Mono vs Naive CD4 T | groups colored after embedding",
  col.grid = "#E6E6E6",
  col.axis = "#555555",
  col.lab = "#333333"
)
graphics::legend(
  "topright",
  legend = c("Up", "Down", "Not significant"),
  col = unname(colors[c("Up", "Down", "Not significant")]),
  pch = 16,
  pt.cex = 1.1,
  bty = "n",
  cex = 0.9
)
grDevices::dev.off()

counts <- table(data$group)
cat("genes_plotted=", nrow(data), "\n", sep = "")
cat("Up=", unname(counts[["Up"]]), "\n", sep = "")
cat("Down=", unname(counts[["Down"]]), "\n", sep = "")
cat("Not_significant=", unname(counts[["Not significant"]]), "\n", sep = "")
cat("pdf=", normalizePath(output_pdf), "\n", sep = "")
cat("umap_3d_check=PASS\n")

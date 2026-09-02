# Plot unsupervised expression-based gene UMAP coordinates calculated in Python.
input_tsv <- "docs/images/pbmc3k-expression-gene-umap.tsv"
output_png <- "docs/images/pbmc3k-expression-gene-umap.png"
output_pdf <- "docs/images/pbmc3k-expression-gene-umap.pdf"

data <- utils::read.delim(input_tsv, check.names = FALSE)
data$group <- factor(
  data$group, levels = c("Up", "Down", "Not significant")
)
colors <- c(
  Up = "#F04438", Down = "#3977B8", `Not significant` = "#B7B7B7"
)

plot <- ggplot2::ggplot(
  data,
  ggplot2::aes(x = UMAP_1, y = UMAP_2, color = group)
) +
  ggplot2::geom_point(size = 0.85, alpha = 0.62, stroke = 0) +
  ggplot2::scale_color_manual(
    values = colors,
    breaks = c("Up", "Down", "Not significant"),
    drop = FALSE,
    name = NULL
  ) +
  ggplot2::labs(
    title = "Expression-based gene UMAP",
    subtitle = "CD14+ Mono vs Naive CD4 T | color added after unsupervised embedding",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_void(base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 11, color = "#666666", hjust = 0.5),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 11),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(18, 24, 18, 24)
  )

stopifnot(capabilities("cairo"))
ggplot2::ggsave(
  output_png, plot = plot, width = 9, height = 6.5,
  units = "in", dpi = 300, type = "cairo"
)
ggplot2::ggsave(
  output_pdf, plot = plot, width = 9, height = 6.5,
  units = "in", device = grDevices::cairo_pdf
)
counts <- table(data$group)
stopifnot(all(is.finite(data$UMAP_1)), all(is.finite(data$UMAP_2)))
cat("genes_plotted=", nrow(data), "\n", sep = "")
cat("Up=", unname(counts[["Up"]]), "\n", sep = "")
cat("Down=", unname(counts[["Down"]]), "\n", sep = "")
cat("Not_significant=", unname(counts[["Not significant"]]), "\n", sep = "")
cat("png=", normalizePath(output_png), "\n", sep = "")
cat("pdf=", normalizePath(output_pdf), "\n", sep = "")
cat("expression_gene_umap_check=PASS\n")

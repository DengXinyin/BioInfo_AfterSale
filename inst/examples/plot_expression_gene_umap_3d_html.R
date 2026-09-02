# Interactive Plotly HTML for rotating the three-dimensional gene UMAP.
if (!requireNamespace("plotly", quietly = TRUE)) stop("plotly is required.")
if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("htmlwidgets is required.")

input_tsv <- "docs/images/pbmc3k-expression-gene-umap-3d.tsv"
output_html <- "docs/images/pbmc3k-expression-gene-umap-3d.html"
data <- utils::read.delim(input_tsv, check.names = FALSE)

colors <- c(Up = "#F04438", Down = "#3977B8", `Not significant` = "#B7B7B7")
sizes <- c(Up = 3.2, Down = 3.2, `Not significant` = 1.8)
opacities <- c(Up = 0.86, Down = 0.88, `Not significant` = 0.26)
plot <- plotly::plot_ly()

# Draw grey points first so highlighted genes remain visible.
for (group in c("Not significant", "Down", "Up")) {
  subset_data <- data[data$group == group, , drop = FALSE]
  hover <- sprintf(
    paste0(
      "Gene: %s<br>Group: %s<br>log2FC: %.3f<br>padj: %.3g",
      "<br>UMAP: %.2f, %.2f, %.2f"
    ),
    subset_data$gene, group, subset_data$avg_log2FC, subset_data$p_val_adj,
    subset_data$UMAP_1, subset_data$UMAP_2, subset_data$UMAP_3
  )
  plot <- plotly::add_trace(
    plot,
    data = subset_data,
    x = ~UMAP_1,
    y = ~UMAP_2,
    z = ~UMAP_3,
    type = "scatter3d",
    mode = "markers",
    name = group,
    text = hover,
    hoverinfo = "text",
    marker = list(
      size = unname(sizes[[group]]),
      color = unname(colors[[group]]),
      opacity = unname(opacities[[group]]),
      line = list(width = 0)
    )
  )
}

plot <- plotly::layout(
  plot,
  title = list(
    text = paste0(
      "3D expression-based gene UMAP",
      "<br><sup>CD14+ Mono vs Naive CD4 T | drag to rotate, scroll to zoom</sup>"
    ),
    x = 0.5
  ),
  legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.05),
  scene = list(
    xaxis = list(title = "UMAP 1", showbackground = TRUE, backgroundcolor = "white"),
    yaxis = list(title = "UMAP 2", showbackground = TRUE, backgroundcolor = "white"),
    zaxis = list(title = "UMAP 3", showbackground = TRUE, backgroundcolor = "white"),
    aspectmode = "data"
  ),
  paper_bgcolor = "white",
  plot_bgcolor = "white",
  margin = list(l = 0, r = 0, b = 40, t = 80)
)

quarto_pandoc <- "/Applications/quarto/bin/tools/aarch64/pandoc"
if (file.exists(quarto_pandoc)) {
  Sys.setenv(RSTUDIO_PANDOC = dirname(quarto_pandoc))
}
htmlwidgets::saveWidget(
  plot,
  file = output_html,
  selfcontained = TRUE,
  title = "3D expression-based gene UMAP"
)
stopifnot(file.exists(output_html), file.info(output_html)$size > 0)
cat("genes_plotted=", nrow(data), "\n", sep = "")
cat("html=", normalizePath(output_html), "\n", sep = "")
cat("interactive_3d_check=PASS\n")

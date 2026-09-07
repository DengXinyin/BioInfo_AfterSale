# Community composition plots consolidated from bar_plot.r, heatmap.r,
# tree_bar_plot.r, fungi_bar.R, fungi_heatmap.R, func2_barplot.R,
# func2_heatmap.R, plot_bar.r, and plot_heatmap.r in micro_v2:v5.42.

micro_abundance_bar <- function(abundance, group = NULL, top_n = 20,
                                relative = TRUE, group_mean = FALSE,
                                colors = NULL, title = "Community composition",
                                font_family = "sans") {
  x <- .micro_matrix(abundance)
  if (relative) {
    totals <- colSums(x)
    if (any(totals <= 0)) stop("All samples must have a positive total.", call. = FALSE)
    x <- sweep(x, 2, totals, "/")
  }
  top_n <- min(as.integer(top_n), nrow(x))
  keep <- names(sort(rowMeans(x), decreasing = TRUE))[seq_len(top_n)]
  other <- colSums(x[setdiff(rownames(x), keep), , drop = FALSE])
  x <- x[keep, , drop = FALSE]
  if (length(setdiff(rownames(abundance), keep))) x <- rbind(x, Others = other)

  groups <- .micro_group(group, colnames(x))
  if (group_mean && !is.null(group)) {
    x <- vapply(levels(groups), function(z) rowMeans(x[, groups == z, drop = FALSE]),
                numeric(nrow(x)))
    if (is.null(dim(x))) x <- matrix(x, ncol = 1, dimnames = list(keep, levels(groups)))
  }
  dat <- .micro_long_matrix(x, "Taxon", "Sample")
  dat$Taxon <- factor(dat$Taxon, levels = rev(rownames(x)))
  if (is.null(colors)) colors <- .micro_palette(nrow(x))

  ggplot2::ggplot(dat, ggplot2::aes(x = Sample, y = Value, fill = Taxon)) +
    ggplot2::geom_col(width = 0.82) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(
      labels = if (relative) function(z) paste0(round(100 * z), "%") else function(z) format(z, big.mark = ","),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = NULL, y = if (relative) "Relative abundance" else "Abundance",
                  title = title) +
    .micro_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

micro_abundance_heatmap <- function(abundance, group = NULL, top_n = 30,
                                    scale = c("row", "none"), cluster_rows = TRUE,
                                    cluster_samples = TRUE,
                                    low = "#2166AC", mid = "white", high = "#B2182B",
                                    title = "Abundance heatmap", font_family = "sans") {
  scale <- match.arg(scale)
  x <- .micro_matrix(abundance)
  keep <- names(sort(rowMeans(x), decreasing = TRUE))[seq_len(min(top_n, nrow(x)))]
  x <- x[keep, , drop = FALSE]
  if (scale == "row") {
    x <- t(scale(t(x)))
    x[!is.finite(x)] <- 0
  }
  if (cluster_rows && nrow(x) > 1) x <- x[stats::hclust(stats::dist(x))$order, , drop = FALSE]
  if (cluster_samples && ncol(x) > 1) x <- x[, stats::hclust(stats::dist(t(x)))$order, drop = FALSE]
  dat <- .micro_long_matrix(x, "Feature", "Sample")
  dat$Feature <- factor(dat$Feature, levels = rev(rownames(x)))
  dat$Sample <- factor(dat$Sample, levels = colnames(x))
  if (!is.null(group)) {
    groups <- .micro_group(group, colnames(abundance))
    dat$Group <- groups[match(dat$Sample, names(groups))]
  }

  ggplot2::ggplot(dat, ggplot2::aes(x = Sample, y = Feature, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0,
                                  name = if (scale == "row") "Z-score" else "Abundance") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    .micro_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.border = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text())
}

micro_clustered_abundance_bar <- function(abundance, ...) {
  x <- .micro_matrix(abundance)
  if (ncol(x) > 1) {
    ord <- stats::hclust(stats::dist(t(sweep(x, 2, pmax(colSums(x), 1), "/"))))$order
    x <- x[, ord, drop = FALSE]
  }
  micro_abundance_bar(x, ..., title = "Cluster-ordered community composition")
}

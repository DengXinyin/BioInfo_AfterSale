# Taxonomic/functional composition, abundance heatmaps and sample correlation.

metage_abundance_bar <- function(abundance, top_n = 20, relative = TRUE,
                                  group = NULL, group_mean = FALSE, colors = NULL,
                                  title = "Taxonomic / functional composition",
                                  font_family = metage_font_family()) {
  x <- .metage_matrix(abundance)
  if (relative) x <- sweep(x, 2, pmax(colSums(x), 1), "/")
  keep <- names(sort(rowMeans(x), decreasing = TRUE))[seq_len(min(top_n, nrow(x)))]
  other_names <- setdiff(rownames(x), keep)
  y <- x[keep, , drop = FALSE]
  if (length(other_names)) y <- rbind(y, Others = colSums(x[other_names, , drop = FALSE]))
  groups <- .metage_group(group, colnames(y))
  if (group_mean) {
    feature_names <- rownames(y)
    y <- do.call(cbind, lapply(levels(groups), function(g) rowMeans(y[, groups == g, drop = FALSE])))
    dimnames(y) <- list(feature_names, levels(groups))
  }
  dat <- .metage_long(y, "Feature", "Sample")
  dat$Feature <- factor(dat$Feature, levels = rev(rownames(y)))
  if (is.null(colors)) colors <- .metage_palette(nrow(y))
  ggplot2::ggplot(dat, ggplot2::aes(Sample, Value, fill = Feature)) +
    ggplot2::geom_col(width = 0.82) + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = if (relative) scales::percent else ggplot2::waiver(), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = if (relative) "Relative abundance" else "Abundance", title = title) +
    .metage_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

metage_abundance_heatmap <- function(abundance, top_n = 30, scale_rows = TRUE,
                                      cluster_rows = TRUE, cluster_samples = TRUE,
                                      title = "Abundance heatmap",
                                      font_family = metage_font_family()) {
  x <- .metage_matrix(abundance)
  keep <- names(sort(rowMeans(x), decreasing = TRUE))[seq_len(min(top_n, nrow(x)))]
  x <- x[keep, , drop = FALSE]
  if (scale_rows) { x <- t(scale(t(x))); x[!is.finite(x)] <- 0 }
  if (cluster_rows && nrow(x) > 1) x <- x[stats::hclust(stats::dist(x))$order, , drop = FALSE]
  if (cluster_samples && ncol(x) > 1) x <- x[, stats::hclust(stats::dist(t(x)))$order, drop = FALSE]
  dat <- .metage_long(x)
  dat$Feature <- factor(dat$Feature, levels = rev(rownames(x)))
  dat$Sample <- factor(dat$Sample, levels = colnames(x))
  ggplot2::ggplot(dat, ggplot2::aes(Sample, Feature, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  name = if (scale_rows) "Z-score" else "Abundance") +
    ggplot2::labs(x = NULL, y = NULL, title = title) + .metage_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.title = ggplot2::element_text(size = 16))
}

metage_sample_correlation <- function(abundance, method = "spearman", cluster = FALSE,
                                       title = "Sample correlation",
                                       font_family = metage_font_family()) {
  x <- .metage_matrix(abundance)
  z <- stats::cor(x, method = method)
  if (cluster && ncol(z) > 1) { o <- stats::hclust(stats::as.dist(1 - z))$order; z <- z[o, o] }
  dat <- .metage_long(z, "Sample1", "Sample2")
  dat$Sample1 <- factor(dat$Sample1, levels = rev(rownames(z)))
  dat$Sample2 <- factor(dat$Sample2, levels = colnames(z))
  ggplot2::ggplot(dat, ggplot2::aes(Sample2, Sample1, fill = Value)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-1, 1), name = "rho") +
    ggplot2::labs(x = NULL, y = NULL, title = title) + .metage_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

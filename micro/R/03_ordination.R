# Ordination and distance visualizations consolidated from PCA.r, PCOA.r,
# NMDS.r, pls_da.r, UPGMA.r, tree.r, and dis_heatmap.r in micro_v2:v5.42.

micro_ordination <- function(abundance, group, method = c("PCoA", "NMDS", "PCA"),
                             distance = "bray", labels = FALSE, ellipse = FALSE,
                             colors = NULL, title = NULL, font_family = "sans",
                             seed = 123) {
  method <- match.arg(method)
  x <- .micro_matrix(abundance)
  groups <- .micro_group(group, colnames(x))
  set.seed(seed)
  subtitle <- NULL
  explained <- c(NA_real_, NA_real_)

  if (method == "PCA") {
    fit <- stats::prcomp(t(x), center = TRUE, scale. = TRUE)
    scores <- as.data.frame(fit$x[, 1:2, drop = FALSE])
    explained <- 100 * fit$sdev[1:2]^2 / sum(fit$sdev^2)
  } else {
    if (!requireNamespace("vegan", quietly = TRUE)) stop("Package `vegan` is required.")
    d <- vegan::vegdist(t(x), method = distance)
    if (method == "PCoA") {
      fit <- stats::cmdscale(d, k = 2, eig = TRUE, add = TRUE)
      scores <- as.data.frame(fit$points)
      names(scores) <- c("Axis1", "Axis2")
      positive <- fit$eig[fit$eig > 0]
      explained <- 100 * fit$eig[1:2] / sum(positive)
    } else {
      fit <- vegan::metaMDS(d, k = 2, trymax = 50, trace = FALSE)
      scores <- as.data.frame(fit$points)
      names(scores) <- c("Axis1", "Axis2")
      subtitle <- paste0("Stress = ", format(round(fit$stress, 3), nsmall = 3))
    }
  }
  names(scores)[1:2] <- c("Axis1", "Axis2")
  scores$Sample <- rownames(scores)
  scores$Group <- groups[match(scores$Sample, names(groups))]
  if (is.null(colors)) colors <- .micro_palette(nlevels(scores$Group))
  xlab <- if (is.finite(explained[1])) sprintf("%s1 (%.1f%%)", method, explained[1]) else paste0(method, "1")
  ylab <- if (is.finite(explained[2])) sprintf("%s2 (%.1f%%)", method, explained[2]) else paste0(method, "2")

  p <- ggplot2::ggplot(scores, ggplot2::aes(x = Axis1, y = Axis2, colour = Group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::labs(x = xlab, y = ylab, title = title %||% method, subtitle = subtitle) +
    .micro_theme(font_family)
  if (labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 3,
                                        max.overlaps = 50, show.legend = FALSE)
    } else p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), size = 3, vjust = -0.6)
  }
  if (ellipse) {
    if (requireNamespace("ggforce", quietly = TRUE)) {
      p <- p + ggforce::geom_mark_ellipse(ggplot2::aes(fill = Group), alpha = 0.08,
                                          show.legend = FALSE)
    } else if (all(table(scores$Group) >= 3)) {
      p <- p + ggplot2::stat_ellipse()
    }
  }
  structure(list(plot = p, scores = scores, model = fit, explained = explained),
            class = "micro_ordination")
}

micro_distance_heatmap <- function(distance, group = NULL, cluster = TRUE,
                                   low = "white", high = "#B2182B",
                                   title = "Sample distance", font_family = "sans") {
  d <- if (inherits(distance, "dist")) as.matrix(distance) else .micro_matrix(distance, "distance")
  if (nrow(d) != ncol(d)) stop("`distance` must be square.")
  if (cluster && nrow(d) > 1) {
    ord <- stats::hclust(stats::as.dist(d))$order
    d <- d[ord, ord, drop = FALSE]
  }
  dat <- .micro_long_matrix(d, "Sample1", "Sample2")
  dat$Sample1 <- factor(dat$Sample1, levels = rev(rownames(d)))
  dat$Sample2 <- factor(dat$Sample2, levels = colnames(d))
  ggplot2::ggplot(dat, ggplot2::aes(x = Sample2, y = Sample1, fill = Value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = low, high = high, name = "Distance") +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL, title = title) +
    .micro_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.border = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text())
}

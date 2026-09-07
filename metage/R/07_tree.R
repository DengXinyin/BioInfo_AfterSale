# A reusable replacement for bar_tree.R; accepts an existing phylogenetic tree.

metage_tree_abundance <- function(tree, abundance, top_n = 30, color = "#377EB8",
                                   font_family = metage_font_family()) {
  .metage_require("ggtree")
  .metage_require("ape")
  if (!inherits(tree, "phylo")) stop("`tree` must be an ape::phylo object.")
  if (is.matrix(abundance) || is.data.frame(abundance)) {
    x <- .metage_matrix(abundance)
    abundance <- rowMeans(x)
  }
  if (!is.numeric(abundance) || is.null(names(abundance))) {
    stop("`abundance` must be a named numeric vector or feature-by-sample matrix.")
  }
  keep <- intersect(names(sort(abundance, decreasing = TRUE)), tree$tip.label)
  keep <- utils::head(keep, top_n)
  if (length(keep) < 2) stop("Fewer than two tree tips match abundance names.")
  tr <- if (length(keep) < length(tree$tip.label)) ape::keep.tip(tree, keep) else tree
  tip <- data.frame(Feature = keep, Abundance = abundance[keep])
  tip$Feature <- factor(tip$Feature, levels = rev(tr$tip.label))
  p_tree <- ggtree::ggtree(tr) + ggtree::geom_tiplab(size = 3, family = font_family)
  p_bar <- ggplot2::ggplot(tip, ggplot2::aes(Abundance, Feature)) +
    ggplot2::geom_col(fill = color, width = 0.7) +
    ggplot2::labs(x = "Mean abundance", y = NULL) + .metage_theme(font_family) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  structure(list(tree = p_tree, abundance = p_bar, data = tip), class = "metage_tree_abundance")
}

metage_tree_abundance_grob <- function(x, widths = c(1.4, 1)) {
  if (!inherits(x, "metage_tree_abundance")) stop("Use `metage_tree_abundance()` first.")
  .metage_require("gridExtra")
  gridExtra::arrangeGrob(x$tree, x$abundance, ncol = 2, widths = widths)
}

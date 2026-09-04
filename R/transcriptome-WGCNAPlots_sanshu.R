#' Draw a WGCNA sample-clustering tree
#'
#' @param expression Numeric sample-by-gene matrix, or an `hclust` object.
#' @param method Linkage method passed to [stats::hclust()].
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional PDF filename.
#' @param width,height Optional output dimensions.
#' @return Invisibly returns the `hclust` object.
#' @export
wgcna_sample_tree <- function(expression, method = "average",
                              title = "Sample clustering",
                              style = NULL, output_file = NULL,
                              width = NULL, height = NULL) {
  style <- .plot_style_or_default(style)
  width <- width %||% style$global$figure_width
  height <- height %||% style$global$figure_height
  .assert_positive_number(width, "width")
  .assert_positive_number(height, "height")
  tree <- if (inherits(expression, "hclust")) {
    expression
  } else {
    matrix <- as.matrix(expression)
    suppressWarnings(storage.mode(matrix) <- "double")
    if (nrow(matrix) < 2L || ncol(matrix) < 2L || any(!is.finite(matrix))) {
      stop("`expression` must contain at least two samples and genes with finite values.", call. = FALSE)
    }
    if (is.null(rownames(matrix))) stop("Expression rows must have sample names.", call. = FALSE)
    stats::hclust(stats::dist(matrix), method = method)
  }
  if (!is.null(output_file)) {
    if (!grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
      stop("`wgcna_sample_tree` writes PDF output only.", call. = FALSE)
    }
    grDevices::cairo_pdf(output_file, width = width, height = height,
                         family = style$global$font_family)
  }
  old <- graphics::par(family = style$global$font_family)
  on.exit({
    graphics::par(old)
    if (!is.null(output_file)) grDevices::dev.off()
  }, add = TRUE)
  graphics::plot(tree, main = title, xlab = "", sub = "",
                 cex.main = style$text$title$size / 18,
                 cex.axis = style$text$axis_text$size / 18)
  invisible(tree)
}

#' Draw WGCNA module gene counts
#'
#' @param data Module count table.
#' @param module,count Column names.
#' @param colors Optional named module colors.
#' @inheritParams wgcna_sample_tree
#' @param dpi Output resolution passed to the shared saver.
#' @return A ggplot object.
#' @export
wgcna_module_sizes <- function(data, module = "Module", count = "Count",
                               colors = NULL, title = "Module gene counts",
                               style = NULL, output_file = NULL, width = NULL,
                               height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(module, count))
  .sanshu_numeric(data, count)
  if (any(data[[count]] < 0) || anyDuplicated(data[[module]])) {
    stop("Modules must be unique and counts non-negative.", call. = FALSE)
  }
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data[order(data[[count]], decreasing = TRUE), , drop = FALSE]
  d[[module]] <- factor(as.character(d[[module]]), levels = as.character(d[[module]]))
  module_levels <- as.character(d[[module]])
  if (is.null(colors)) {
    candidate <- tolower(module_levels)
    valid <- !inherits(try(grDevices::col2rgb(candidate), silent = TRUE), "try-error")
    colors <- if (valid) stats::setNames(candidate, module_levels) else NULL
  }
  palette <- .sanshu_palette(module_levels, colors, z$style)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[module]], y = .data[[count]], fill = .data[[module]]
  )) +
    ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = .data[[count]]), vjust = -0.3,
                       family = z$style$global$font_family,
                       size = z$style$text$data_label$size / 3.2) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::labs(x = "Module", y = "Gene count", title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      family = z$style$global$font_family, angle = 45, hjust = 1
    ))
  .sanshu_finish(p, z)
}

#' Draw a WGCNA module-trait correlation heatmap
#'
#' @param correlation,pvalue Matrices with identical dimensions and names.
#' @param title Plot title.
#' @inheritParams wgcna_module_sizes
#' @return A ggplot object.
#' @export
wgcna_module_trait_heatmap <- function(correlation, pvalue,
                                       title = "Module-trait relationships",
                                       style = NULL, output_file = NULL,
                                       width = NULL, height = NULL, dpi = NULL) {
  cor_matrix <- as.matrix(correlation)
  p_matrix <- as.matrix(pvalue)
  suppressWarnings(storage.mode(cor_matrix) <- "double")
  suppressWarnings(storage.mode(p_matrix) <- "double")
  if (!identical(dim(cor_matrix), dim(p_matrix)) ||
      !identical(dimnames(cor_matrix), dimnames(p_matrix))) {
    stop("Correlation and P-value matrices must have identical dimensions and names.", call. = FALSE)
  }
  if (!length(cor_matrix) || any(!is.finite(cor_matrix)) || any(!is.finite(p_matrix)) ||
      any(cor_matrix < -1 | cor_matrix > 1) || any(p_matrix < 0 | p_matrix > 1)) {
    stop("Matrices contain invalid correlation or P values.", call. = FALSE)
  }
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE)
  names(d) <- c("Module", "Trait", "Correlation")
  d$Pvalue <- as.vector(p_matrix)
  d$Module <- factor(d$Module, levels = rev(rownames(cor_matrix)))
  d$Trait <- factor(d$Trait, levels = colnames(cor_matrix))
  d$Label <- sprintf("%.2f\n(p=%.2g)", d$Correlation, d$Pvalue)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$Trait, .data$Module,
                                       fill = .data$Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label),
                       family = z$style$global$font_family,
                       size = z$style$text$data_label$size / 3.2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  limits = c(-1, 1), midpoint = 0,
                                  breaks = c(-1, 0, 1), name = "Correlation") +
    ggplot2::labs(x = "Trait", y = "Module", title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      family = z$style$global$font_family, angle = 45, hjust = 1
    ), panel.grid = ggplot2::element_blank())
  .sanshu_finish(p, z)
}

#' Draw WGCNA module-membership versus gene-significance scatter
#'
#' @param data Gene-level WGCNA statistics.
#' @param module_membership,gene_significance Column names.
#' @param module Optional module column used for color.
#' @param colors Optional module colors.
#' @param title Plot title.
#' @inheritParams wgcna_module_sizes
#' @return A ggplot object.
#' @export
wgcna_mm_gs <- function(data, module_membership = "ModuleMembership",
                        gene_significance = "GeneSignificance", module = NULL,
                        colors = NULL, title = "Module membership vs gene significance",
                        style = NULL, output_file = NULL, width = NULL,
                        height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(module_membership, gene_significance, module))
  .sanshu_numeric(data, c(module_membership, gene_significance))
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  d$.Module <- if (is.null(module)) "Genes" else as.character(d[[module]])
  levels <- unique(d$.Module)
  palette <- .sanshu_palette(levels, colors, z$style)
  correlation <- stats::cor(d[[module_membership]], d[[gene_significance]], method = "pearson")
  pvalue <- stats::cor.test(d[[module_membership]], d[[gene_significance]])$p.value
  subtitle <- sprintf("Pearson r = %.3f; P = %.3g", correlation, pvalue)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[module_membership]], y = .data[[gene_significance]], color = .data$.Module
  )) +
    ggplot2::geom_point(alpha = 0.65, size = 1.5) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         color = "black", fill = "#BDBDBD", linewidth = 0.7) +
    ggplot2::scale_color_manual(values = palette, breaks = levels,
                                name = if (is.null(module)) NULL else "Module",
                                guide = if (is.null(module)) "none" else "legend") +
    ggplot2::labs(x = "Module membership", y = "Gene significance",
                  title = title, subtitle = subtitle) +
    z$style$ggplot_theme
  .sanshu_finish(p, z)
}

.plot_style_or_default <- function(style) {
  if (is.null(style)) style <- choose_plot_style()
  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by choose_plot_style().", call. = FALSE)
  }
  style
}

.save_bioinfo_plot <- function(plot, style, output_file, width, height, dpi) {
  if (is.null(output_file)) return(invisible(plot))
  ggplot2::ggsave(
    output_file, plot,
    width = width %||% style$global$figure_width,
    height = height %||% style$global$figure_height,
    dpi = dpi %||% style$global$dpi, bg = "white"
  )
  invisible(plot)
}

#' Plot PCA sample scores
#'
#' @param data Data frame containing sample scores and groups.
#' @param pc1,pc2,group,sample Column names for PC1, PC2, group, and labels.
#' @param variance Optional two-element percentage vector for axis labels.
#' @param ellipse Draw group confidence ellipses.
#' @param ellipse_level Confidence level.
#' @param show_labels Display sample labels.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
plot_pca <- function(data, pc1 = "PC1", pc2 = "PC2", group = "Group",
                     sample = "Sample", variance = NULL, ellipse = TRUE,
                     ellipse_level = 0.95, show_labels = TRUE,
                     title = "Principal Component Analysis", style = NULL,
                     output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  style <- .plot_style_or_default(style)
  required <- c(pc1, pc2, group, if (show_labels) sample)
  if (length(miss <- setdiff(required, names(data))))
    stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d <- as.data.frame(data); d[[group]] <- factor(d[[group]])
  pal <- stats::setNames(rep(style$group_palette, length.out = nlevels(d[[group]])), levels(d[[group]]))
  xlab <- if (is.null(variance)) pc1 else sprintf("%s (%.2f%%)", pc1, variance[1])
  ylab <- if (is.null(variance)) pc2 else sprintf("%s (%.2f%%)", pc2, variance[2])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[pc1]], y = .data[[pc2]], color = .data[[group]]))
  if (ellipse) p <- p + ggplot2::stat_ellipse(ggplot2::aes(fill = .data[[group]]), geom = "polygon", level = ellipse_level, alpha = .16, linewidth = .7)
  p <- p + ggplot2::geom_point(size = 3) + ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_fill_manual(values = pal) + ggplot2::labs(x = xlab, y = ylab, title = title, color = NULL, fill = NULL)
  if (show_labels) p <- p + ggplot2::geom_text(ggplot2::aes(label = .data[[sample]]), vjust = -0.7, check_overlap = TRUE, show.legend = FALSE)
  p <- p + style$ggplot_theme
  .save_bioinfo_plot(p, style, output_file, width, height, dpi); p
}

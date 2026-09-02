#' Plot a circular cross-correlation heatmap
#' @inheritParams plot_dot_heatmap
#' @return A ggplot object using polar coordinates.
#' @export
plot_circular_heatmap <- function(x, y, method = "spearman", title = NULL, style = NULL, output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  z <- .validate_style_output(style, output_file, width, height, dpi); rho <- .correlation_matrix(x, y, method); d <- .correlation_long(rho); d$row <- factor(d$row, levels = rev(rownames(rho))); d$column <- factor(d$column, levels = colnames(rho))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$column, .data$row, fill = .data$value)) + ggplot2::geom_tile(color = "white", linewidth = .25, na.rm = FALSE) + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-1, 1), na.value = "grey90", name = paste0(method, " r")) + ggplot2::coord_polar(clip = "off") + ggplot2::labs(x = NULL, y = NULL, title = title) + z$style$ggplot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(size = z$style$text$axis_text$size * .65), axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  .save_bioinfo_plot(p, z$style, output_file, z$width, z$height, z$dpi); p
}

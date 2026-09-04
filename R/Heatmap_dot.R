#' Plot a dot correlation heatmap
#' @param x,y Numeric data frames with samples in rows and variables in columns.
#' @param method Correlation method passed to [stats::cor()].
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object with the correlation matrix in `plot$data`.
#' @export
plot_dot_heatmap <- function(x, y, method = "spearman", title = NULL, style = NULL, output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  z <- .validate_style_output(style, output_file, width, height, dpi); d <- .correlation_long(.correlation_matrix(x, y, method))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$column, .data$row)) + ggplot2::geom_point(ggplot2::aes(size = .data$size, color = .data$value), na.rm = TRUE) + ggplot2::scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-1, 1), name = paste0(method, " r")) + ggplot2::scale_size_continuous(range = c(1, 10), limits = c(0, 1), name = paste0("|", method, " r|")) + ggplot2::labs(x = NULL, y = NULL, title = title) + z$style$ggplot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1), panel.grid = ggplot2::element_line(color = "#E1E1E1", linewidth = .3))
  .save_bioinfo_plot(p, z$style, output_file, z$width, z$height, z$dpi); p
}

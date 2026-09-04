#' Plot a rectangular cross-correlation heatmap
#' @inheritParams plot_dot_heatmap
#' @return A ggplot object with a tile-based correlation matrix.
#' @export
plot_cross_heatmap <- function(x, y, method = "spearman", title = NULL, style = NULL, output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  z <- .validate_style_output(style, output_file, width, height, dpi); stat <- .correlation_stats(x, y, method); rho <- stat$rho; d <- .correlation_long(rho); d$stars <- .significance_stars(as.vector(stat$p))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$column, .data$row, fill = .data$value)) + ggplot2::geom_tile(color = "white", linewidth = .2) + ggplot2::geom_text(ggplot2::aes(label = .data$stars), size = z$style$text$data_label$size / 3.2, na.rm = TRUE) + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-1, 1), na.value = "grey90", name = paste0(method, " r")) + ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = "* p < 0.05; ** p < 0.01; *** p < 0.001") + z$style$ggplot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1), panel.grid = ggplot2::element_blank(), aspect.ratio = max(0.35, min(1.6, nrow(rho) / ncol(rho))))
  .save_bioinfo_plot(p, z$style, output_file, z$width, z$height, z$dpi); p
}

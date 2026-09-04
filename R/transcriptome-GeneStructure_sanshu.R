#' Draw differential exon usage with a gene model
#'
#' @param data Long-form exon expression table.
#' @param exon,group,expression,padj Column names.
#' @param exon_order Explicit exon order; defaults to first appearance.
#' @param colors Optional group colors.
#' @param highlight_color Color for the most significant exon.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A gtable containing expression and gene-model panels.
#' @export
deu_exon_expression <- function(
    data, exon = "Exon", group = "Group", expression = "Expression",
    padj = "Padj", exon_order = NULL, colors = NULL,
    highlight_color = "#FF00FF", title = "Differential exon usage",
    style = NULL, output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(exon, group, expression, padj))
  .sanshu_numeric(data, c(expression, padj))
  if (any(data[[padj]] < 0 | data[[padj]] > 1)) stop("Padj must be between 0 and 1.", call. = FALSE)
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for the two-panel exon plot.", call. = FALSE)
  }
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  observed <- unique(as.character(d[[exon]]))
  if (is.null(exon_order)) exon_order <- observed
  if (!setequal(exon_order, observed) || anyDuplicated(exon_order)) {
    stop("`exon_order` must contain every exon exactly once.", call. = FALSE)
  }
  d[[exon]] <- factor(as.character(d[[exon]]), levels = exon_order)
  groups <- unique(as.character(d[[group]]))
  palette <- .sanshu_palette(groups, colors, z$style)
  top <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[exon]], y = .data[[expression]], color = .data[[group]],
    group = .data[[group]]
  )) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = palette, breaks = groups, name = "Group") +
    ggplot2::labs(x = NULL, y = "Expression", title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  exon_padj <- tapply(d[[padj]], d[[exon]], min)
  highlight <- names(which.min(exon_padj))
  model <- data.frame(
    Exon = factor(exon_order, levels = exon_order), Y = 1,
    Type = ifelse(exon_order == highlight, "Most significant DEU exon", "Other exon")
  )
  model_colors <- c("Most significant DEU exon" = highlight_color, "Other exon" = "#999999")
  bottom <- ggplot2::ggplot(model, ggplot2::aes(.data$Exon, .data$Y)) +
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "black", linewidth = 0.6) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Type), width = 0.75, height = 0.45,
                       color = "black", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = model_colors, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.5)) +
    ggplot2::labs(x = "Exon", y = NULL) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank())
  arranged <- gridExtra::arrangeGrob(top, bottom, ncol = 1, heights = c(3, 1))
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, arranged, width = z$width, height = z$height,
                    dpi = z$dpi, device = grDevices::cairo_pdf)
  }
  arranged
}

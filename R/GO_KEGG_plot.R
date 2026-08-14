#' Plot an enrichment result
#'
#' Filters an enrichment result and creates a dot plot or bar plot with
#' `enrichplot`. Common font, title, legend, and axis-text settings are exposed
#' for routine after-sales figure changes. The returned object is a regular
#' ggplot and can be modified further with ggplot2.
#'
#' @param result An `enrichResult` or `gseaResult` object.
#' @param plot_type Either `"dotplot"` or `"barplot"`.
#' @param filter_by Column used to filter terms: `"p.adjust"` or `"pvalue"`.
#' @param cutoff Maximum value retained from `filter_by`.
#' @param show_category Maximum number of categories displayed.
#' @param x X variable accepted by [enrichplot::dotplot()], normally
#'   `"GeneRatio"` or `"Count"`. Ignored for bar plots.
#' @param color Color variable passed to enrichplot.
#' @param order_by Ordering variable passed to enrichplot.
#' @param decreasing Whether ordering is decreasing.
#' @param title Optional plot title.
#' @param font_family Font family.
#' @param base_size Base text size.
#' @param title_size Optional title size.
#' @param axis_text_size Optional axis tick-label size.
#' @param legend_text_size Optional legend text size.
#' @param legend_position A ggplot2 legend position.
#' @param label_format Maximum label width used by enrichplot.
#' @param ... Additional arguments passed to the enrichplot function.
#'
#' @return A ggplot object.
#' @export
GO_KEGG_plot <- function(
    result,
    plot_type = c("dotplot", "barplot"),
    filter_by = c("p.adjust", "pvalue"),
    cutoff = 0.05,
    show_category = 20,
    x = "GeneRatio",
    color = filter_by,
    order_by = x,
    decreasing = TRUE,
    title = NULL,
    font_family = "Arial",
    base_size = 14,
    title_size = base_size + 2,
    axis_text_size = base_size,
    legend_text_size = max(base_size - 2, 1),
    legend_position = "right",
    label_format = 40,
    ...) {

  plot_type <- match.arg(plot_type)
  filter_by <- match.arg(filter_by)
  .assert_probability(cutoff, "cutoff")
  .assert_positive_number(show_category, "show_category")
  .assert_positive_number(base_size, "base_size")

  filtered <- .filter_result_object(result, filter_by, cutoff)
  if (!nrow(as.data.frame(filtered))) {
    stop(sprintf("No enrichment terms pass `%s <= %s`.", filter_by, cutoff), call. = FALSE)
  }

  if (plot_type == "dotplot") {
    plot <- enrichplot::dotplot(
      filtered, x = x, color = color, showCategory = show_category,
      orderBy = order_by, decreasing = decreasing,
      label_format = label_format, ...
    )
  } else {
    plot <- graphics::barplot(
      filtered, color = color, showCategory = show_category, ...
    )
  }

  plot +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family, size = title_size, hjust = 0.5
      ),
      axis.title = ggplot2::element_text(family = font_family, size = base_size),
      axis.text = ggplot2::element_text(
        family = font_family, size = axis_text_size, color = "black"
      ),
      legend.title = ggplot2::element_text(family = font_family, size = base_size),
      legend.text = ggplot2::element_text(
        family = font_family, size = legend_text_size
      ),
      legend.position = legend_position,
      panel.grid.minor = ggplot2::element_blank()
    )
}

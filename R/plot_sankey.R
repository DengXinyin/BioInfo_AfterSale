#' Draw a Sankey diagram
#'
#' Draws an alluvial/Sankey diagram from a long-form edge table (`source`,
#' `target`, `value`) using [ggalluvial]. Columns are drawn as strata and
#' flows sized by the value.
#'
#' @param data Data frame with `source`, `target`, and `value` columns.
#' @param source_column Name of the source (left) column.
#' @param target_column Name of the target (right) column.
#' @param value_column Name of the numeric value column.
#' @param fill_column Column used to colour the strata/flows. `NULL` uses
#'   `target_column`.
#' @param curve_type Strand/curve geometry passed to
#'   [ggalluvial::geom_flow()]: `"linear"`, `"cubic"`, `"quintic"`, `"sine"`,
#'   `"arctangent"`, `"sigmoid"`, or `"xspline"`.
#' @param flow_alpha Flow fill opacity.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @param x_label,y_label Axis labels.
#' @param legend_title Legend title.
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()].
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object.
#' @export
plot_sankey <- function(
    data,
    source_column = "source",
    target_column = "target",
    value_column = "value",
    fill_column = NULL,
    curve_type = c("linear", "sigmoid", "cubic", "quintic", "sine",
                   "arctangent", "xspline"),
    flow_alpha = 0.6,
    title = "Sankey Diagram",
    subtitle = NULL,
    x_label = "Omics Type",
    y_label = "Gene/Protein Count",
    legend_title = "Pathway",
    font_family = "sans",
    style = NULL,
    output_file = NULL,
    figure_width = 10,
    figure_height = 8,
    dpi = 300) {

  curve_type <- match.arg(curve_type)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (column in c(source_column, target_column, value_column)) {
    if (!is.character(column) || length(column) != 1L || is.na(column) ||
        !nzchar(column) || !column %in% names(data)) {
      stop("Column `", column, "` is not a column in `data`.", call. = FALSE)
    }
  }
  if (!is.numeric(data[[value_column]])) {
    stop("`value_column` must be numeric.", call. = FALSE)
  }
  .assert_probability(flow_alpha, "flow_alpha")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  if (is.null(fill_column)) fill_column <- target_column
  if (!is.character(fill_column) || length(fill_column) != 1L ||
      !fill_column %in% names(data)) {
    stop("`fill_column` must be a column in `data`.", call. = FALSE)
  }

  if (is.null(style)) {
    style <- choose_plot_style(
      font_family = font_family,
      title = list(size = 16),
      subtitle = list(size = 13),
      axis_title = list(size = 13),
      axis_text = list(size = 11),
      legend_title = list(size = 12),
      legend_text = list(size = 11),
      legend = list(position = "right")
    )
  }
  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
  }

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[source_column]],
      stratum = .data[[target_column]],
      alluvium = interaction(.data[[source_column]], .data[[target_column]]),
      y = .data[[value_column]],
      fill = .data[[fill_column]],
      label = .data[[target_column]]
    )
  ) +
    ggalluvial::geom_flow(alpha = flow_alpha, curve_type = curve_type) +
    ggalluvial::geom_stratum(alpha = 0.8) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, size = 3.5,
                       family = font_family) +
    ggplot2::scale_fill_viridis_d(option = "D", name = legend_title) +
    ggplot2::labs(x = x_label, y = y_label, title = title,
                  subtitle = subtitle) +
    style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

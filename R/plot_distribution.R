#' Draw a distribution histogram, density, or both
#'
#' Plots the distribution of a numeric column as a histogram, a density curve,
#' or an overlay of both, with optional grouping by a second column. This
#' reimplements the 08_分布图 module of the 组学可视化百宝箱 using only
#' [ggplot2].
#'
#' @param data Data frame containing the value column.
#' @param value_column Name of the numeric value column.
#' @param group_column Optional name of the grouping column used to split the
#'   distribution by fill/colour. `NULL` draws a single overall distribution.
#' @param group_colors Named or unnamed colour vector mapped to the group
#'   levels. `NULL` uses the toolbox palette.
#' @param plot_type Which layers to draw: `"histogram"`, `"density"`, or
#'   `"both"`.
#' @param bins Number of histogram bins.
#' @param fill_color Fill colour used when there is no grouping.
#' @param line_color Line colour used when there is no grouping.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels. `NULL` uses the value column name and a
#'   sensible default y label.
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()].
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object.
#' @export
plot_distribution <- function(
    data,
    value_column = "Expression",
    group_column = NULL,
    group_colors = NULL,
    plot_type = c("both", "histogram", "density"),
    bins = 20,
    fill_color = "#4DBBD5",
    line_color = "#3C8DBC",
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    font_family = "sans",
    style = NULL,
    output_file = NULL,
    figure_width = 10,
    figure_height = 8,
    dpi = 300) {

  plot_type <- match.arg(plot_type)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(value_column) || length(value_column) != 1L ||
      is.na(value_column) || !nzchar(value_column) ||
      !value_column %in% names(data)) {
    stop("Column `", value_column, "` is not a column in `data`.", call. = FALSE)
  }
  if (!is.numeric(data[[value_column]])) {
    stop("`value_column` must be numeric.", call. = FALSE)
  }
  .assert_positive_number(bins, "bins")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  bins <- as.integer(bins)
  if (bins < 1L) bins <- 1L
  if (!is.character(fill_color) || length(fill_color) != 1L ||
      is.na(fill_color) || !nzchar(fill_color)) {
    stop("`fill_color` must be one non-empty colour value.", call. = FALSE)
  }
  if (!is.character(line_color) || length(line_color) != 1L ||
      is.na(line_color) || !nzchar(line_color)) {
    stop("`line_color` must be one non-empty colour value.", call. = FALSE)
  }

  if (is.null(x_label)) x_label <- value_column
  if (is.null(y_label)) y_label <- if (plot_type == "density") "Density" else "Count"

  if (is.null(style)) {
    style <- choose_plot_style(
      font_family = font_family,
      title = list(size = 16),
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

  has_group <- !is.null(group_column)
  if (has_group) {
    if (!is.character(group_column) || length(group_column) != 1L ||
        !group_column %in% names(data)) {
      stop("`group_column` must be a column in `data`.", call. = FALSE)
    }
    data[[group_column]] <- if (is.factor(data[[group_column]])) {
      droplevels(data[[group_column]])
    } else {
      factor(as.character(data[[group_column]]),
             unique(as.character(data[[group_column]])))
    }
    group_levels <- levels(data[[group_column]])
    if (is.null(group_colors)) {
      palette <- .toolbox_group_palette()
      if (length(group_levels) > length(palette)) {
        group_colors <- grDevices::hcl.colors(length(group_levels), "Dark 3")
      } else {
        group_colors <- palette[seq_along(group_levels)]
      }
      names(group_colors) <- group_levels
    }
    if (is.null(names(group_colors)) || anyNA(group_colors) ||
        any(!nzchar(group_colors))) {
      stop("`group_colors` must be a named character vector.", call. = FALSE)
    }
    missing_colors <- setdiff(group_levels, names(group_colors))
    if (length(missing_colors)) {
      stop("`group_colors` is missing level(s): ",
           paste(missing_colors, collapse = ", "), call. = FALSE)
    }
    base_aes <- ggplot2::aes(
      x = .data[[value_column]], fill = .data[[group_column]],
      colour = .data[[group_column]]
    )
  } else {
    base_aes <- ggplot2::aes(x = .data[[value_column]])
  }

  plot <- ggplot2::ggplot(data, base_aes)

  if (plot_type %in% c("histogram", "both")) {
    hist_aes <- ggplot2::aes(y = ggplot2::after_stat(density))
    plot <- plot + ggplot2::geom_histogram(
      mapping = if (plot_type == "both") hist_aes else NULL,
      bins = bins, alpha = if (has_group) 0.5 else 0.8, colour = "white"
    )
  }
  if (plot_type %in% c("density", "both") || has_group) {
    plot <- plot + ggplot2::geom_density(alpha = 0.4, linewidth = 0.9)
  }

  if (has_group) {
    plot <- plot +
      ggplot2::scale_fill_manual(
        values = group_colors, breaks = group_levels, drop = FALSE,
        name = group_column
      ) +
      ggplot2::scale_colour_manual(
        values = group_colors, breaks = group_levels, drop = FALSE,
        name = group_column
      )
  } else {
    plot <- plot +
      ggplot2::scale_fill_manual(
        values = stats::setNames(fill_color, "overall"), guide = "none"
      ) +
      ggplot2::scale_colour_manual(
        values = stats::setNames(line_color, "overall"), guide = "none"
      )
  }

  plot <- plot +
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    style$ggplot_theme

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

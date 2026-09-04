# Toolbox qualitative palette (matches the 组学可视化百宝箱 four/five-colour
# convention used by violin/box plots).
.toolbox_qualitative_palette <- function() {
  c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
  )
}

#' Draw a violin or box plot with optional jitter
#'
#' Plots a long-form value distribution split by a grouping column using
#' [ggplot2::geom_violin()], [ggplot2::geom_boxplot()], and an optional
#' [ggplot2::geom_jitter()] overlay. Supports optional faceting when the data
#' contains a feature/measurement column.
#'
#' @param data Data frame in long format with a value column and a group column.
#' @param value_column Name of the numeric value column.
#' @param group_column Name of the grouping column used on the x axis.
#' @param fill_column Optional name of the column mapped to fill. `NULL` uses
#'   `group_column`.
#' @param facet_column Optional name of a column used by
#'   [ggplot2::facet_wrap()]. `NULL` draws a single panel.
#' @param group_colors Named or unnamed colour vector mapped to the group levels.
#'   `NULL` uses the toolbox qualitative palette.
#' @param plot_type Which layers to draw: `"violin"`, `"box"`, or `"both"`.
#' @param show_jitter Whether to overlay a [ggplot2::geom_jitter()] layer.
#' @param jitter_alpha Jitter point opacity.
#' @param jitter_width Horizontal jitter width.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels. `NULL` uses the group and value column
#'   names.
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()]. When
#'   supplied, it replaces the font and text-size defaults.
#' @param facet_ncol Number of facet columns when `facet_column` is set.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object.
#' @export
plot_violin_box <- function(
    data,
    value_column = "Value",
    group_column = "Species",
    fill_column = NULL,
    facet_column = NULL,
    group_colors = NULL,
    plot_type = c("both", "violin", "box"),
    show_jitter = TRUE,
    jitter_alpha = 0.3,
    jitter_width = 0.08,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    font_family = "Times New Roman",
    style = NULL,
    facet_ncol = 2,
    output_file = NULL,
    figure_width = 10,
    figure_height = 8,
    dpi = 300) {

  plot_type <- match.arg(plot_type)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (column in c(value_column, group_column)) {
    if (!is.character(column) || length(column) != 1L || is.na(column) ||
        !nzchar(column) || !column %in% names(data)) {
      stop("Column `", column, "` is not a column in `data`.", call. = FALSE)
    }
  }
  if (!is.numeric(data[[value_column]])) {
    stop("`value_column` must be numeric.", call. = FALSE)
  }
  if (is.null(fill_column)) fill_column <- group_column
  if (!is.character(fill_column) || length(fill_column) != 1L ||
      !fill_column %in% names(data)) {
    stop("`fill_column` must be a column in `data`.", call. = FALSE)
  }
  .assert_flag(show_jitter, "show_jitter")
  .assert_probability(jitter_alpha, "jitter_alpha")
  .assert_positive_number(jitter_width, "jitter_width")
  .assert_positive_number(facet_ncol, "facet_ncol")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  facet_ncol <- as.integer(facet_ncol)
  if (facet_ncol < 1L) facet_ncol <- 1L

  data[[group_column]] <- if (is.factor(data[[group_column]])) {
    droplevels(data[[group_column]])
  } else {
    factor(as.character(data[[group_column]]), unique(as.character(data[[group_column]])))
  }
  group_levels <- levels(data[[group_column]])
  fill_levels <- if (identical(fill_column, group_column)) {
    group_levels
  } else {
    factor(data[[fill_column]]); unique(as.character(data[[fill_column]]))
  }

  if (is.null(group_colors)) {
    palette <- .toolbox_qualitative_palette()
    if (length(fill_levels) > length(palette)) {
      group_colors <- grDevices::hcl.colors(length(fill_levels), "Dark 3")
    } else {
      group_colors <- palette[seq_along(fill_levels)]
    }
    names(group_colors) <- fill_levels
  }
  if (is.null(names(group_colors)) || anyNA(group_colors) ||
      any(!nzchar(group_colors))) {
    stop("`group_colors` must be a named character vector.", call. = FALSE)
  }
  missing_colors <- setdiff(fill_levels, names(group_colors))
  if (length(missing_colors)) {
    stop("`group_colors` is missing level(s): ",
         paste(missing_colors, collapse = ", "), call. = FALSE)
  }
  if (inherits(try(grDevices::col2rgb(group_colors), silent = TRUE), "try-error")) {
    stop("`group_colors` must contain valid R colors.", call. = FALSE)
  }

  if (is.null(x_label)) x_label <- group_column
  if (is.null(y_label)) y_label <- value_column

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

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data[[group_column]], y = .data[[value_column]],
                 fill = .data[[fill_column]])
  )

  if (plot_type %in% c("violin", "both")) {
    plot <- plot + ggplot2::geom_violin(trim = FALSE, alpha = 0.6)
  }
  if (plot_type %in% c("box", "both")) {
    plot <- plot + ggplot2::geom_boxplot(
      width = 0.12, fill = "white", outlier.shape = NA
    )
  }
  if (show_jitter) {
    plot <- plot + ggplot2::geom_jitter(
      width = jitter_width, alpha = jitter_alpha, size = 1.5
    )
  }

  plot <- plot +
    ggplot2::scale_fill_manual(
      values = group_colors, breaks = fill_levels, drop = FALSE, name = fill_column
    ) +
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    style$ggplot_theme

  if (!is.null(facet_column)) {
    if (!is.character(facet_column) || length(facet_column) != 1L ||
        !facet_column %in% names(data)) {
      stop("`facet_column` must be a column in `data`.", call. = FALSE)
    }
    plot <- plot + ggplot2::facet_wrap(
      stats::as.formula(paste("~", facet_column)),
      scales = "free_y", ncol = facet_ncol
    ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

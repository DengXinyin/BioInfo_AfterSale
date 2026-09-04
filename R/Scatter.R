#' Draw a scatter plot with trend line and correlation annotation
#'
#' Plots two numeric columns against each other with an optional linear or
#' LOESS trend line and an optional per-panel correlation annotation via
#' [ggpubr::stat_cor()]. Supports optional grouping, faceting, and marginal
#' density/box layers via [ggExtra::ggMarginal()].
#'
#' @param data Data frame containing the x and y columns.
#' @param x_column Name of the numeric x column.
#' @param y_column Name of the numeric y column.
#' @param group_column Optional name of the grouping column. `NULL` draws a
#'   single overall panel.
#' @param facet_column Optional name of a column used by
#'   [ggplot2::facet_wrap()]. `NULL` draws a single panel.
#' @param group_colors Named or unnamed colour vector mapped to the group
#'   levels. `NULL` uses the toolbox palette.
#' @param cor_method Correlation coefficient to annotate: `"pearson"` or
#'   `"spearman"` (passed to [ggpubr::stat_cor()]). `NULL` hides the annotation.
#' @param smooth_method Trend-line smoother: `"lm"`, `"loess"`, or `NULL` for
#'   no trend line.
#' @param smooth_se Whether to draw the confidence band.
#' @param point_color Point colour used when there is no grouping.
#' @param point_size Point size.
#' @param point_alpha Point opacity.
#' @param marginal Whether to add marginal density/box layers via
#'   [ggExtra::ggMarginal()]. `"none"`, `"density"`, or `"boxplot"`.
#' @param marginal_type Type passed to [ggExtra::ggMarginal()] when `marginal`
#'   is not `"none"`.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels. `NULL` uses the column names.
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()].
#' @param facet_ncol Number of facet columns when `facet_column` is set.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object (or the [ggExtra::ggMarginal()] result when
#'   `marginal` is not `"none"`).
#' @export
plot_scatter <- function(
    data,
    x_column = "VarX",
    y_column = "VarY",
    group_column = NULL,
    facet_column = NULL,
    group_colors = NULL,
    cor_method = c("pearson", "spearman"),
    smooth_method = c("lm", "loess"),
    smooth_se = TRUE,
    point_color = "#2C3E50",
    point_size = 2,
    point_alpha = 0.6,
    marginal = c("none", "density", "boxplot"),
    marginal_type = "density",
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    font_family = "Times New Roman",
    style = NULL,
    facet_ncol = 3,
    output_file = NULL,
    figure_width = 10,
    figure_height = 8,
    dpi = 300) {

  cor_method <- match.arg(cor_method)
  smooth_method <- match.arg(smooth_method)
  marginal <- match.arg(marginal)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (column in c(x_column, y_column)) {
    if (!is.character(column) || length(column) != 1L || is.na(column) ||
        !nzchar(column) || !column %in% names(data)) {
      stop("Column `", column, "` is not a column in `data`.", call. = FALSE)
    }
    if (!is.numeric(data[[column]])) {
      stop("Column `", column, "` must be numeric.", call. = FALSE)
    }
  }
  .assert_positive_number(point_size, "point_size")
  .assert_probability(point_alpha, "point_alpha")
  .assert_flag(smooth_se, "smooth_se")
  .assert_positive_number(facet_ncol, "facet_ncol")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  if (!is.character(point_color) || length(point_color) != 1L ||
      is.na(point_color) || !nzchar(point_color)) {
    stop("`point_color` must be one non-empty colour value.", call. = FALSE)
  }
  facet_ncol <- as.integer(facet_ncol)
  if (facet_ncol < 1L) facet_ncol <- 1L

  if (is.null(x_label)) x_label <- x_column
  if (is.null(y_label)) y_label <- y_column

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

  if (!is.null(group_column)) {
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
      x = .data[[x_column]], y = .data[[y_column]], color = .data[[group_column]]
    )
  } else {
    base_aes <- ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]])
  }

  plot <- ggplot2::ggplot(data, base_aes) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha)

  if (!is.null(group_column)) {
    plot <- plot +
      ggplot2::scale_color_manual(
        values = group_colors, breaks = group_levels, drop = FALSE,
        name = group_column
      ) +
      ggplot2::scale_fill_manual(
        values = group_colors, breaks = group_levels, drop = FALSE,
        name = group_column
      )
  } else {
    plot <- plot +
      ggplot2::scale_color_manual(values = stats::setNames(point_color, "overall"))
  }

  if (!is.null(smooth_method)) {
    smooth_args <- list(method = smooth_method, se = smooth_se, alpha = 0.15)
    if (!is.null(group_column)) {
      plot <- plot + ggplot2::geom_smooth(
        method = smooth_method, se = smooth_se, alpha = 0.15,
        ggplot2::aes(fill = .data[[group_column]])
      )
    } else {
      smooth_args$color <- point_color
      smooth_args$fill <- grDevices::adjustcolor(point_color, alpha.f = 0.3)
      plot <- plot + do.call(ggplot2::geom_smooth, smooth_args)
    }
  }

  if (!is.null(cor_method)) {
    plot <- plot + ggpubr::stat_cor(
      method = cor_method, size = 3.5, color = "black"
    )
  }

  plot <- plot +
    ggplot2::labs(x = x_label, y = y_label, title = title, colour = group_column) +
    style$ggplot_theme
  if (!is.null(group_column)) {
    plot <- plot + ggplot2::theme(legend.position = "bottom")
  }

  if (!is.null(facet_column)) {
    if (!is.character(facet_column) || length(facet_column) != 1L ||
        !facet_column %in% names(data)) {
      stop("`facet_column` must be a column in `data`.", call. = FALSE)
    }
    plot <- plot + ggplot2::facet_wrap(
      stats::as.formula(paste("~", facet_column)), ncol = facet_ncol
    )
  }

  if (marginal != "none") {
    plot <- ggExtra::ggMarginal(
      plot, type = marginal_type, groupColour = !is.null(group_column),
      groupFill = !is.null(group_column), alpha = 0.3
    )
  }

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

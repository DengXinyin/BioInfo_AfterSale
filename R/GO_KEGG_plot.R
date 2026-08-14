#' Plot an enrichment result
#'
#' Filters an enrichment result and creates a dot plot or bar plot. Native
#' `enrichResult` objects use `enrichplot`; standardized data frames use a
#' flexible ggplot2 path that supports arbitrary numeric axes and colors.
#'
#' @param result An `enrichResult`, `gseaResult`, or standardized data frame.
#' @param plot_type Either `"dotplot"` or `"barplot"`.
#' @param filter_by Column used to filter terms: `"p.adjust"`, `"pvalue"`, or
#'   `NULL` to keep every supplied row.
#' @param cutoff Maximum value retained from `filter_by`.
#' @param show_category Maximum number of categories displayed.
#' @param x X variable accepted by [enrichplot::dotplot()], normally
#'   `"GeneRatio"` or `"Count"`. For data frames, any numeric column.
#' @param x_transform X-axis transformation for data frames: `"identity"` or
#'   `"neg_log10"`.
#' @param color Color variable passed to enrichplot.
#' @param label Label column used for data-frame plots.
#' @param size Point-size column used for data-frame dot plots.
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
#' @param x_label,y_label Optional axis labels. Labels can be character values
#'   or plotmath expressions such as `expression(-log[10](Pvalue))`.
#' @param color_label,size_label Optional legend labels for data-frame plots.
#' @param color_palette Continuous colors used for data-frame plots.
#' @param size_breaks Optional numeric Count legend breaks for data-frame dot
#'   plots. `NULL` uses ggplot2 automatic breaks.
#' @param size_range Two positive numbers defining the minimum and maximum
#'   point sizes.
#' @param point_alpha Point opacity from 0 to 1.
#' @param x_limits,x_breaks Optional numeric X-axis limits and breaks.
#' @param x_expand Optional numeric expansion passed to
#'   [ggplot2::scale_x_continuous()].
#' @param legend_order Character vector defining guide order. Supported values
#'   are `"color"` and `"size"`.
#' @param style Optional object returned by [choose_plot_style()]. When
#'   supplied, it replaces the individual legacy font and legend arguments.
#' @param output_file Optional output filename. When supplied, the function
#'   saves the plot and still returns the ggplot object.
#' @param figure_width,figure_height Output width and height in inches. `NULL`
#'   uses the values stored in `style`.
#' @param dpi Output resolution. `NULL` uses the value stored in `style`.
#' @param device Optional graphics device passed to [ggplot2::ggsave()]. PDF
#'   files default to [grDevices::cairo_pdf()].
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
    x_transform = c("identity", "neg_log10"),
    color = filter_by,
    label = "Description",
    size = "Count",
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
    x_label = NULL,
    y_label = NULL,
    color_label = NULL,
    size_label = NULL,
    color_palette = c("#2166AC", "#00BFC4", "#7AD151", "#FDE725", "#D73027"),
    size_breaks = NULL,
    size_range = c(4, 12),
    point_alpha = 0.9,
    x_limits = NULL,
    x_breaks = NULL,
    x_expand = NULL,
    legend_order = c("color", "size"),
    style = NULL,
    output_file = NULL,
    figure_width = NULL,
    figure_height = NULL,
    dpi = NULL,
    device = NULL,
    ...) {

  color_was_missing <- missing(color)
  plot_type <- match.arg(plot_type)
  x_transform <- match.arg(x_transform)
  if (!is.null(filter_by)) filter_by <- match.arg(filter_by)
  if (color_was_missing && is.null(filter_by)) color <- x
  .assert_probability(cutoff, "cutoff")
  .assert_probability(point_alpha, "point_alpha")
  .assert_positive_number(show_category, "show_category")
  .assert_positive_number(base_size, "base_size")
  .assert_numeric_vector(size_range, "size_range", length_required = 2L, positive = TRUE)
  if (size_range[[1]] > size_range[[2]]) {
    stop("`size_range` must be ordered from minimum to maximum.", call. = FALSE)
  }
  if (!is.null(size_breaks)) {
    .assert_numeric_vector(size_breaks, "size_breaks", positive = TRUE)
  }
  if (!is.null(x_limits)) {
    .assert_numeric_vector(x_limits, "x_limits", length_required = 2L)
    if (x_limits[[1]] >= x_limits[[2]]) {
      stop("`x_limits` must be ordered from minimum to maximum.", call. = FALSE)
    }
  }
  if (!is.null(x_breaks)) .assert_numeric_vector(x_breaks, "x_breaks")
  if (!is.null(x_expand)) .assert_numeric_vector(x_expand, "x_expand")
  if (!is.character(legend_order) || !length(legend_order) ||
      anyNA(legend_order) || anyDuplicated(legend_order) ||
      any(!legend_order %in% c("color", "size"))) {
    stop(
      "`legend_order` must contain unique values selected from 'color' and 'size'.",
      call. = FALSE
    )
  }
  guide_order <- stats::setNames(seq_along(legend_order), legend_order)
  color_guide_order <- if ("color" %in% names(guide_order)) {
    unname(guide_order[["color"]])
  } else 0
  size_guide_order <- if ("size" %in% names(guide_order)) {
    unname(guide_order[["size"]])
  } else 0
  if (is.null(style)) {
    style <- choose_plot_style(
      font_family = font_family,
      title = list(size = title_size),
      axis_title = list(size = base_size),
      axis_text = list(size = axis_text_size),
      legend_title = list(size = base_size),
      legend_text = list(size = legend_text_size),
      legend = list(
        show = !identical(legend_position, "none"),
        position = legend_position
      )
    )
  }
  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
  }

  if (is.data.frame(result)) {
    table <- as.data.frame(result, stringsAsFactors = FALSE)
    required <- unique(c(label, x, color, if (plot_type == "dotplot") size))
    if (!is.null(filter_by)) required <- unique(c(required, filter_by))
    missing_columns <- setdiff(required, names(table))
    if (length(missing_columns)) {
      stop(
        "Data-frame result is missing column(s): ",
        paste(missing_columns, collapse = ", "), call. = FALSE
      )
    }
    if (!is.null(filter_by)) {
      filter_values <- suppressWarnings(as.numeric(table[[filter_by]]))
      table <- table[!is.na(filter_values) & filter_values <= cutoff, , drop = FALSE]
    }
    if (!nrow(table)) {
      stop("No enrichment terms remain for plotting.", call. = FALSE)
    }

    x_values <- suppressWarnings(as.numeric(table[[x]]))
    if (x_transform == "neg_log10") {
      x_values <- -log10(pmax(x_values, .Machine$double.xmin))
    }
    color_values <- suppressWarnings(as.numeric(table[[color]]))
    size_values <- if (plot_type == "dotplot") {
      suppressWarnings(as.numeric(table[[size]]))
    } else {
      rep(NA_real_, nrow(table))
    }
    keep <- is.finite(x_values) & is.finite(color_values)
    if (plot_type == "dotplot") keep <- keep & is.finite(size_values)
    table <- table[keep, , drop = FALSE]
    x_values <- x_values[keep]
    color_values <- color_values[keep]
    size_values <- size_values[keep]
    if (!nrow(table)) stop("Plot mappings contain no finite numeric values.", call. = FALSE)

    labels <- vapply(as.character(table[[label]]), function(value) {
      paste(strwrap(value, width = label_format), collapse = "\n")
    }, character(1))
    order_values <- if (identical(order_by, x)) {
      x_values
    } else if (order_by %in% names(table)) {
      suppressWarnings(as.numeric(table[[order_by]]))
    } else {
      stop("`order_by` is not present in the data-frame result.", call. = FALSE)
    }
    ordering <- order(order_values, decreasing = decreasing, na.last = NA)
    if (length(ordering) > show_category) ordering <- ordering[seq_len(show_category)]
    table <- table[ordering, , drop = FALSE]
    table$.XValue <- x_values[ordering]
    table$.ColorValue <- color_values[ordering]
    table$.SizeValue <- size_values[ordering]
    table$.Label <- factor(labels[ordering], levels = rev(unique(labels[ordering])))

    if (is.null(x_label)) {
      x_label <- if (x_transform == "neg_log10") {
        expression(-log[10](Pvalue))
      } else {
        x
      }
    }
    if (is.null(color_label)) color_label <- color
    if (is.null(size_label)) size_label <- size

    if (plot_type == "dotplot") {
      size_scale_args <- list(name = size_label, range = size_range)
      if (!is.null(size_breaks)) size_scale_args$breaks <- size_breaks
      plot <- ggplot2::ggplot(
        table,
        ggplot2::aes(x = .XValue, y = .Label, color = .ColorValue, size = .SizeValue)
      ) +
        ggplot2::geom_point(alpha = point_alpha) +
        do.call(ggplot2::scale_size_continuous, size_scale_args) +
        ggplot2::scale_color_gradientn(colors = color_palette, name = color_label) +
        ggplot2::guides(
          color = ggplot2::guide_colorbar(order = color_guide_order),
          size = ggplot2::guide_legend(order = size_guide_order)
        )
    } else {
      plot <- ggplot2::ggplot(
        table, ggplot2::aes(x = .XValue, y = .Label, fill = .ColorValue)
      ) +
        ggplot2::geom_col(width = 0.75) +
        ggplot2::scale_fill_gradientn(colors = color_palette, name = color_label) +
        ggplot2::guides(
          fill = ggplot2::guide_colorbar(order = color_guide_order)
        )
    }
    plot <- plot +
      ggplot2::labs(x = x_label, y = y_label, title = title) +
      style$ggplot_theme
    x_scale_args <- list()
    if (!is.null(x_limits)) x_scale_args$limits <- x_limits
    if (!is.null(x_breaks)) x_scale_args$breaks <- x_breaks
    if (!is.null(x_expand)) x_scale_args$expand <- x_expand
    if (length(x_scale_args)) {
      plot <- plot + do.call(ggplot2::scale_x_continuous, x_scale_args)
    }
    return(.finish_enrichment_plot(
      plot, style, output_file, figure_width, figure_height, dpi, device
    ))
  }

  if (!methods::is(result, "enrichResult") && !methods::is(result, "gseaResult")) {
    stop("`result` must be an enrichResult, gseaResult, or data.frame.", call. = FALSE)
  }
  filtered <- if (is.null(filter_by)) {
    result
  } else {
    .filter_result_object(result, filter_by, cutoff)
  }
  if (!nrow(as.data.frame(filtered))) stop("No enrichment terms remain for plotting.", call. = FALSE)

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

  plot <- plot + ggplot2::labs(title = title) + style$ggplot_theme
  if (!is.null(x_label)) plot <- plot + ggplot2::labs(x = x_label)
  if (!is.null(y_label)) plot <- plot + ggplot2::labs(y = y_label)
  x_scale_args <- list()
  if (!is.null(x_limits)) x_scale_args$limits <- x_limits
  if (!is.null(x_breaks)) x_scale_args$breaks <- x_breaks
  if (!is.null(x_expand)) x_scale_args$expand <- x_expand
  if (length(x_scale_args)) {
    plot <- plot + do.call(ggplot2::scale_x_continuous, x_scale_args)
  }
  .finish_enrichment_plot(
    plot, style, output_file, figure_width, figure_height, dpi, device
  )
}

.finish_enrichment_plot <- function(
    plot, style, output_file, figure_width, figure_height, dpi, device) {
  figure_width <- figure_width %||% style$global$figure_width
  figure_height <- figure_height %||% style$global$figure_height
  dpi <- dpi %||% style$global$dpi
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1L ||
        is.na(output_file) || !nzchar(output_file)) {
      stop("`output_file` must be one non-empty filename.", call. = FALSE)
    }
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      stop("Output directory does not exist: ", output_dir, call. = FALSE)
    }
    if (is.null(device) && identical(tolower(tools::file_ext(output_file)), "pdf")) {
      device <- grDevices::cairo_pdf
    }
    save_args <- list(
      filename = output_file,
      plot = plot,
      width = figure_width,
      height = figure_height,
      units = "in",
      dpi = dpi
    )
    if (!is.null(device)) save_args$device <- device
    do.call(ggplot2::ggsave, save_args)
  }
  plot
}

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
#' @param cutoff Strict upper cutoff for `filter_by`; retained rows satisfy
#'   `0 < filter_by < cutoff`. This excludes zero-valued P values because they
#'   cannot be represented faithfully by `-log10()`.
#' @param show_category Maximum number of categories displayed.
#' @param x X variable accepted by [enrichplot::dotplot()], normally
#'   `"GeneRatio"` or `"Count"`. For data frames, any numeric column.
#' @param x_transform X-axis transformation for data frames: `"identity"` or
#'   `"neg_log10"`.
#' @param color Color variable passed to enrichplot.
#' @param color_transform Color transformation for data-frame plots:
#'   `"identity"` or `"neg_log10"`.
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
#' @param color_breaks Optional numeric vector of at most three color-legend
#'   breaks. By default, three breaks are selected from the displayed values.
#' @param color_palette Continuous colors used for data-frame plots.
#' @param size_breaks Optional numeric Count legend breaks for data-frame dot
#'   plots. `NULL` selects up to three representative observed values.
#' @param size_range Two positive numbers defining the minimum and maximum
#'   point sizes.
#' @param point_alpha Point opacity from 0 to 1.
#' @param highlight_terms Optional character vector of pathway or term names to
#'   emphasize on the Y axis. Values are matched exactly against `label`.
#' @param highlight_color Color used for highlighted term labels.
#' @param highlight_bold Whether highlighted term labels are bold.
#' @param x_limits,x_breaks Optional numeric X-axis limits and breaks.
#'   `x_limits` filters the plotted data; use `x_zoom` for visual zooming
#'   without removing observations.
#' @param x_zoom Optional numeric X-axis zoom passed to
#'   [ggplot2::coord_cartesian()].
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
# All GO/KEGG plot text prefers Times New Roman.
GO_KEGG_plot <- function(
    result,
    plot_type = c("dotplot", "barplot"),
    filter_by = c("p.adjust", "pvalue"),
    cutoff = 0.05,
    show_category = 20,
    x = "GeneRatio",
    x_transform = c("identity", "neg_log10"),
    color = filter_by,
    color_transform = c("identity", "neg_log10"),
    label = "Description",
    size = "Count",
    order_by = x,
    decreasing = TRUE,
    title = NULL,
    font_family = "Times New Roman",
    base_size = 14,
    title_size = base_size + 2,
    axis_text_size = base_size,
    legend_text_size = max(base_size - 2, 1),
    legend_position = "right",
    label_format = 40,
    x_label = NULL,
    y_label = NULL,
    color_label = NULL,
    color_breaks = NULL,
    size_label = NULL,
    color_palette = c("#2166AC", "#00BFC4", "#7AD151", "#FDE725", "#D73027"),
    size_breaks = NULL,
    size_range = c(4, 12),
    point_alpha = 0.9,
    highlight_terms = NULL,
    highlight_color = "#D73027",
    highlight_bold = TRUE,
    x_limits = NULL,
    x_breaks = NULL,
    x_expand = NULL,
    x_zoom = NULL,
    legend_order = c("color", "size"),
    style = NULL,
    output_file = NULL,
    figure_width = NULL,
    figure_height = NULL,
    dpi = NULL,
    device = NULL,
    ...) {

  # Enrichment/differential-result plotting uses the strict rule 0 < p < 0.05
  # (or the supplied cutoff). A reported p = 0 is treated as underflow or a
  # below-detection value, excluded from filtering, and never plotted as an
  # exact finite probability.

  color_was_missing <- missing(color)
  plot_type <- match.arg(plot_type)
  x_transform <- match.arg(x_transform)
  color_transform <- match.arg(color_transform)
  if (!is.null(filter_by)) filter_by <- match.arg(filter_by)
  if (color_was_missing && is.null(filter_by)) {
    color <- if (is.data.frame(result)) x else "p.adjust"
  }
  .assert_probability(cutoff, "cutoff")
  .assert_probability(point_alpha, "point_alpha")
  .assert_flag(highlight_bold, "highlight_bold")
  .assert_positive_number(show_category, "show_category")
  .assert_positive_number(base_size, "base_size")
  .assert_numeric_vector(size_range, "size_range", length_required = 2L, positive = TRUE)
  if (size_range[[1]] > size_range[[2]]) {
    stop("`size_range` must be ordered from minimum to maximum.", call. = FALSE)
  }
  if (!is.null(size_breaks)) {
    .assert_numeric_vector(size_breaks, "size_breaks", positive = TRUE)
  }
  if (!is.null(color_breaks)) {
    .assert_numeric_vector(color_breaks, "color_breaks", positive = FALSE)
    if (length(color_breaks) > 3L) stop("`color_breaks` may contain at most three values.", call. = FALSE)
  }
  if (!is.null(highlight_terms)) {
    if (!is.character(highlight_terms) || anyNA(highlight_terms) ||
        any(!nzchar(highlight_terms))) {
      stop("`highlight_terms` must be NULL or a non-empty character vector.", call. = FALSE)
    }
    highlight_terms <- unique(highlight_terms)
  }
  if (!is.character(highlight_color) || length(highlight_color) != 1L ||
      is.na(highlight_color) || !nzchar(highlight_color)) {
    stop("`highlight_color` must be one non-empty color value.", call. = FALSE)
  }
  highlight_rgb <- tryCatch(
    grDevices::col2rgb(highlight_color),
    error = function(error) NULL
  )
  if (is.null(highlight_rgb)) {
    stop("`highlight_color` must be a valid R color.", call. = FALSE)
  }
  highlight_color <- grDevices::rgb(
    highlight_rgb[1, 1], highlight_rgb[2, 1], highlight_rgb[3, 1],
    maxColorValue = 255
  )
  if (!is.null(x_limits)) {
    .assert_numeric_vector(x_limits, "x_limits", length_required = 2L)
    if (x_limits[[1]] >= x_limits[[2]]) {
      stop("`x_limits` must be ordered from minimum to maximum.", call. = FALSE)
    }
  }
  if (!is.null(x_breaks)) .assert_numeric_vector(x_breaks, "x_breaks")
  if (!is.null(x_expand)) .assert_numeric_vector(x_expand, "x_expand")
  if (!is.null(x_zoom)) {
    .assert_numeric_vector(x_zoom, "x_zoom", length_required = 2L)
    if (x_zoom[[1]] >= x_zoom[[2]]) {
      stop("`x_zoom` must be ordered from minimum to maximum.", call. = FALSE)
    }
  }
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
      keep <- .filter_positive_pvalues(table[[filter_by]], filter_by, cutoff)
      table <- table[keep, , drop = FALSE]
    }
    if (!nrow(table)) {
      stop("No enrichment terms remain for plotting.", call. = FALSE)
    }

    x_values <- suppressWarnings(as.numeric(table[[x]]))
    if (x_transform == "neg_log10") {
      non_positive <- !is.na(x_values) & x_values <= 0
      if (any(non_positive)) {
        warning(
          sprintf(
            "%d row(s) with %s <= 0 were excluded before -log10 transformation.",
            sum(non_positive), x
          ),
          call. = FALSE
        )
      }
      x_values[non_positive] <- NA_real_
      x_values <- -log10(x_values)
    }
    color_values <- suppressWarnings(as.numeric(table[[color]]))
    if (color_transform == "neg_log10") {
      color_values[color_values <= 0] <- NA_real_
      color_values <- -log10(color_values)
    }
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

    term_names <- as.character(table[[label]])
    labels <- vapply(term_names, function(value) {
      paste(strwrap(value, width = label_format), collapse = "\n")
    }, character(1))
    order_values <- if (identical(order_by, x) ||
                        (identical(x_transform, "neg_log10") &&
                         order_by %in% c("pvalue", "p.adjust"))) {
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
    display_labels <- labels[ordering]
    displayed_terms <- term_names[ordering]
    use_highlight_labels <- !is.null(highlight_terms)
    if (use_highlight_labels) {
      highlighted <- displayed_terms %in% highlight_terms
      missing_highlights <- setdiff(highlight_terms, displayed_terms)
      if (length(missing_highlights)) {
        warning(
          "Highlighted term(s) are not displayed: ",
          paste(missing_highlights, collapse = ", "), call. = FALSE
        )
      }
      display_labels <- gsub(
        "\n", "<br>", .escape_html(display_labels), fixed = TRUE
      )
      emphasis <- paste0(
        "color:", highlight_color, ";",
        if (highlight_bold) "font-weight:bold;" else ""
      )
      display_labels[highlighted] <- paste0(
        "<span style=\"", emphasis, "\">",
        display_labels[highlighted], "</span>"
      )
    }
    table$.Label <- factor(
      display_labels, levels = rev(unique(display_labels))
    )

    if (is.null(x_label)) {
      x_label <- if (x_transform == "neg_log10") {
        if (identical(x, "p.adjust")) expression(-log[10](Padj)) else expression(-log[10](Pvalue))
      } else {
        x
      }
    }
    if (is.null(color_label)) {
      color_label <- if (color_transform == "neg_log10" && identical(color, "p.adjust")) {
        expression(-log[10](Padj))
      } else if (color_transform == "neg_log10" && identical(color, "pvalue")) {
        expression(-log[10](Pvalue))
      } else color
    }
    if (is.null(color_breaks)) color_breaks <- .default_numeric_breaks(table$.ColorValue, n = 3L)
    if (is.null(size_label)) size_label <- size

    if (plot_type == "dotplot") {
      size_scale_args <- list(name = size_label, range = size_range)
      size_scale_args$breaks <- if (is.null(size_breaks)) {
        .default_size_breaks(table$.SizeValue, n = 3L)
      } else {
        size_breaks
      }
      plot <- ggplot2::ggplot(
        table,
        ggplot2::aes(x = .XValue, y = .Label, color = .ColorValue, size = .SizeValue)
      ) +
        ggplot2::geom_point(alpha = point_alpha) +
        do.call(ggplot2::scale_size_continuous, size_scale_args) +
        ggplot2::scale_color_gradientn(colors = color_palette, name = color_label, breaks = color_breaks) +
        ggplot2::guides(
          color = ggplot2::guide_colorbar(order = color_guide_order),
          size = ggplot2::guide_legend(order = size_guide_order)
        )
    } else {
      plot <- ggplot2::ggplot(
        table, ggplot2::aes(x = .XValue, y = .Label, fill = .ColorValue)
      ) +
        ggplot2::geom_col(width = 0.75) +
        ggplot2::scale_fill_gradientn(colors = color_palette, name = color_label, breaks = color_breaks) +
        ggplot2::guides(
          fill = ggplot2::guide_colorbar(order = color_guide_order)
        )
    }
    plot <- plot +
      ggplot2::labs(x = x_label, y = y_label, title = title) +
      style$ggplot_theme
    if (use_highlight_labels && style$text$axis_text$show) {
      axis_text_style <- style$text$axis_text
      axis_text_family <- if (nzchar(axis_text_style$font_family)) {
        axis_text_style$font_family
      } else {
        style$global$font_family
      }
      plot <- plot + ggplot2::theme(
        axis.text.y = ggtext::element_markdown(
          family = axis_text_family,
          size = axis_text_style$size,
          face = .text_face(axis_text_style),
          hjust = c(left = 0, center = 0.5, right = 1)[[axis_text_style$align]],
          color = "black"
        )
      )
    }
    x_scale_args <- list()
    if (!is.null(x_limits)) x_scale_args$limits <- x_limits
    if (!is.null(x_breaks)) x_scale_args$breaks <- x_breaks
    if (!is.null(x_expand)) x_scale_args$expand <- x_expand
    if (length(x_scale_args)) {
      plot <- plot + do.call(ggplot2::scale_x_continuous, x_scale_args)
    }
    if (!is.null(x_zoom)) plot <- plot + ggplot2::coord_cartesian(xlim = x_zoom)
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
    native_size_values <- if (size %in% names(plot$data)) {
      suppressWarnings(as.numeric(plot$data[[size]]))
    } else {
      numeric()
    }
    native_size_breaks <- if (is.null(size_breaks)) {
      .default_size_breaks(native_size_values, n = 3L)
    } else {
      size_breaks
    }
    if (length(native_size_breaks)) {
      plot <- suppressMessages(
        plot + ggplot2::scale_size_continuous(
          name = size_label %||% size,
          range = size_range,
          breaks = native_size_breaks
        )
      )
    }
    native_color_values <- if (color %in% names(plot$data)) {
      suppressWarnings(as.numeric(plot$data[[color]]))
    } else numeric()
    native_color_label <- color_label %||% if (identical(color, "p.adjust")) {
      expression(-log[10](Padj))
    } else if (identical(color, "pvalue")) {
      expression(-log[10](Pvalue))
    } else color
    if (color %in% c("p.adjust", "pvalue") && length(native_color_values)) {
      native_color_values[native_color_values <= 0] <- NA_real_
      plot$data$.ColorValue <- -log10(native_color_values)
      plot$mapping$colour <- ggplot2::aes(.data$.ColorValue)$colour
      native_color_values <- plot$data$.ColorValue
    }
    if (length(native_color_values)) {
      native_color_breaks <- color_breaks %||% .default_numeric_breaks(native_color_values, n = 3L)
      plot <- plot + ggplot2::scale_color_gradientn(
        colors = color_palette, name = native_color_label,
        breaks = native_color_breaks
      )
    }
  } else {
    plot <- graphics::barplot(
      filtered, color = color, showCategory = show_category, ...
    )
  }

  plot <- plot + ggplot2::labs(title = title) + style$ggplot_theme
  if (!is.null(highlight_terms) && plot_type == "dotplot") {
    displayed_terms <- if ("Description" %in% names(plot$data)) {
      unique(as.character(plot$data$Description))
    } else {
      character()
    }
    missing_highlights <- setdiff(highlight_terms, displayed_terms)
    if (length(missing_highlights)) {
      warning(
        "Highlighted term(s) are not displayed: ",
        paste(missing_highlights, collapse = ", "), call. = FALSE
      )
    }
    label_formatter <- function(values) {
      wrapped <- vapply(values, function(value) {
        paste(strwrap(value, width = label_format), collapse = "\n")
      }, character(1))
      formatted <- gsub("\n", "<br>", .escape_html(wrapped), fixed = TRUE)
      highlighted <- values %in% highlight_terms
      emphasis <- paste0(
        "color:", highlight_color, ";",
        if (highlight_bold) "font-weight:bold;" else ""
      )
      formatted[highlighted] <- paste0(
        "<span style=\"", emphasis, "\">",
        formatted[highlighted], "</span>"
      )
      formatted
    }
    plot <- suppressMessages(
      plot + ggplot2::scale_y_discrete(labels = label_formatter)
    )
    if (style$text$axis_text$show) {
      axis_text_style <- style$text$axis_text
      axis_text_family <- if (nzchar(axis_text_style$font_family)) {
        axis_text_style$font_family
      } else {
        style$global$font_family
      }
      plot <- plot + ggplot2::theme(
        axis.text.y = ggtext::element_markdown(
          family = axis_text_family,
          size = axis_text_style$size,
          face = .text_face(axis_text_style),
          hjust = c(left = 0, center = 0.5, right = 1)[[axis_text_style$align]],
          color = "black"
        )
      )
    }
  }
  if (!is.null(x_label)) plot <- plot + ggplot2::labs(x = x_label)
  if (!is.null(y_label)) plot <- plot + ggplot2::labs(y = y_label)
  x_scale_args <- list()
  if (!is.null(x_limits)) x_scale_args$limits <- x_limits
  if (!is.null(x_breaks)) x_scale_args$breaks <- x_breaks
  if (!is.null(x_expand)) x_scale_args$expand <- x_expand
  if (length(x_scale_args)) {
    plot <- plot + do.call(ggplot2::scale_x_continuous, x_scale_args)
  }
  if (!is.null(x_zoom)) plot <- plot + ggplot2::coord_cartesian(xlim = x_zoom)
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

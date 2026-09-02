#' Draw a differential-expression volcano plot
#'
#' Creates a three-class volcano plot from a differential-expression result
#' table. Significant genes satisfy `padj > 0`, `padj < padj_cutoff`, and
#' `abs(log2FoldChange) >= log2fc_cutoff`.
#'
#' Rows with missing, non-finite, or non-positive adjusted P values are removed
#' before the `-log10(padj)` transformation. In particular, `padj = 0` is
#' treated as a numerical underflow/below-detection value rather than plotted
#' as an exact probability.
#'
#' @param result Data frame containing differential-expression results.
#' @param log2fc_column Column containing log2 fold changes.
#' @param padj_column Column containing adjusted P values.
#' @param log2fc_cutoff Absolute log2 fold-change threshold. Defaults to 1.
#' @param padj_cutoff Adjusted P-value threshold. Defaults to 0.05.
#' @param group_colors Named colors for `Up`, `Down`, and `Not significant`.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels. Defaults to plotmath expressions for
#'   log2 fold change and -log10 adjusted P value.
#' @param font_family Font family used by the shared plot style. `"sans"` is
#'   the portable default.
#' @param style Optional object returned by [choose_plot_style()]. When
#'   supplied, it replaces the font and text-size defaults.
#' @param point_size Point size for genes.
#' @param point_alpha Point opacity from 0 to 1.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object with the classified data stored in `plot$data`.
#' @export
volcano_plot <- function(
    result,
    log2fc_column = "log2FoldChange",
    padj_column = "padj",
    log2fc_cutoff = 1,
    padj_cutoff = 0.05,
    group_colors = c(
      Up = "#D73027",
      Down = "#4575B4",
      `Not significant` = "#9E9E9E"
    ),
    title = NULL,
    x_label = expression(log[2](FoldChange)),
    y_label = expression(-log[10](padj)),
    font_family = "sans",
    style = NULL,
    point_size = 1.8,
    point_alpha = 0.75,
    output_file = NULL,
    figure_width = 8,
    figure_height = 7,
    dpi = 300) {

  if (!is.data.frame(result)) {
    stop("`result` must be a data frame.", call. = FALSE)
  }
  columns <- c(log2fc_column, padj_column)
  if (!is.character(columns) || anyNA(columns) || any(!nzchar(columns))) {
    stop("`log2fc_column` and `padj_column` must be non-empty names.", call. = FALSE)
  }
  missing_columns <- setdiff(columns, names(result))
  if (length(missing_columns)) {
    stop(
      "Result is missing column(s): ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  .assert_positive_number(log2fc_cutoff, "log2fc_cutoff")
  .assert_probability(padj_cutoff, "padj_cutoff")
  if (padj_cutoff <= 0) stop("`padj_cutoff` must be greater than 0.", call. = FALSE)
  .assert_positive_number(point_size, "point_size")
  .assert_probability(point_alpha, "point_alpha")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")

  required_groups <- c("Up", "Down", "Not significant")
  if (!is.character(group_colors) || is.null(names(group_colors)) ||
      anyNA(group_colors) || any(!nzchar(group_colors)) ||
      any(!required_groups %in% names(group_colors))) {
    stop(
      "`group_colors` must be a named color vector containing Up, Down, and Not significant.",
      call. = FALSE
    )
  }
  group_colors <- group_colors[required_groups]
  if (inherits(try(grDevices::col2rgb(group_colors), silent = TRUE), "try-error")) {
    stop("`group_colors` must contain valid R colors.", call. = FALSE)
  }

  data <- data.frame(
    .log2FoldChange = suppressWarnings(as.numeric(result[[log2fc_column]])),
    .padj = suppressWarnings(as.numeric(result[[padj_column]])),
    stringsAsFactors = FALSE
  )
  valid <- is.finite(data$.log2FoldChange) & is.finite(data$.padj) & data$.padj > 0
  invalid_count <- sum(!valid)
  if (invalid_count) {
    warning(
      sprintf(
        "%d row(s) with missing, non-finite, or non-positive %s were excluded.",
        invalid_count, padj_column
      ),
      call. = FALSE
    )
  }
  data <- data[valid, , drop = FALSE]
  if (!nrow(data)) stop("No valid rows remain for the volcano plot.", call. = FALSE)
  data$.minus_log10_padj <- -log10(data$.padj)
  data$.Group <- "Not significant"
  significant <- data$.padj < padj_cutoff
  data$.Group[significant & data$.log2FoldChange >= log2fc_cutoff] <- "Up"
  data$.Group[significant & data$.log2FoldChange <= -log2fc_cutoff] <- "Down"
  data$.Group <- factor(data$.Group, levels = required_groups)

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
    ggplot2::aes(x = .log2FoldChange, y = .minus_log10_padj, color = .Group)
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_vline(
      xintercept = c(-log2fc_cutoff, log2fc_cutoff),
      linetype = "dashed", color = "black", linewidth = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed", color = "black", linewidth = 0.45
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      breaks = required_groups,
      drop = FALSE,
      name = "Group"
    ) +
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    style$ggplot_theme

  attr(plot, "log2fc_cutoff") <- log2fc_cutoff
  attr(plot, "padj_cutoff") <- padj_cutoff
  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .volcano_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

.volcano_save_plot <- function(plot, output_file, width, height, dpi) {
  if (!is.character(output_file) || length(output_file) != 1L ||
      is.na(output_file) || !nzchar(output_file)) {
    stop("`output_file` must be one non-empty filename.", call. = FALSE)
  }
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir, call. = FALSE)
  }
  extension <- tolower(tools::file_ext(output_file))
  if (!extension %in% c("pdf", "png")) {
    stop("`output_file` must end in .pdf or .png.", call. = FALSE)
  }
  save_args <- list(
    filename = output_file, plot = plot, width = width, height = height,
    units = "in", dpi = dpi, limitsize = FALSE
  )
  if (extension == "pdf") {
    save_args$device <- grDevices::cairo_pdf
  } else {
    save_args$device <- function(filename, width, height, ...) {
      grDevices::png(
        filename = filename, width = width, height = height,
        units = "in", res = dpi, type = "cairo"
      )
    }
  }
  do.call(ggplot2::ggsave, save_args)
  invisible(NULL)
}

#' Draw sequencing error-rate profiles
#'
#' @param data Long-form read-position data.
#' @param position,value,stage,sample Column names for read position, error
#'   rate, filtering stage, and sample ID. `value` is expected as a percentage.
#' @param read Optional read-direction column. When supplied, positions for the
#'   second read are shifted after the first read within each panel.
#' @param colors Named stage colors.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
read_qc_error <- function(data, position = "Position", value = "ErrorRate",
                          stage = "Stage", sample = "SampleID", read = NULL,
                          colors = NULL, title = "Sequencing error rate",
                          style = NULL, output_file = NULL, width = NULL,
                          height = NULL, dpi = NULL) {
  columns <- c(position, value, stage, sample, read)
  columns <- columns[!is.null(columns)]
  .sanshu_require_columns(data, columns)
  .sanshu_numeric(data, c(position, value))
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  if (!is.null(read)) {
    key <- interaction(d[[sample]], d[[stage]], drop = TRUE)
    first <- !grepl("2$|read2|r2", tolower(as.character(d[[read]])))
    maxima <- tapply(d[[position]][first], key[first], max)
    shift <- unname(maxima[as.character(key)])
    second <- !first
    if (any(!is.finite(shift[second]))) {
      stop("Every second-read panel must have matching first-read positions.", call. = FALSE)
    }
    d[[position]][second] <- d[[position]][second] + shift[second]
  }
  stages <- unique(as.character(d[[stage]]))
  palette <- .sanshu_palette(stages, colors, z$style)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[position]], y = .data[[value]], fill = .data[[stage]]
  )) +
    ggplot2::geom_col(width = 0.9, show.legend = FALSE) +
    ggplot2::facet_grid(stats::as.formula(paste(sample, "~", stage)), scales = "free_x") +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::labs(x = "Position", y = "Error rate (%)", title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  .sanshu_finish(p, z)
}

#' Draw base-content profiles before and after filtering
#'
#' @inheritParams read_qc_error
#' @param base Column containing A, T, C, G, N, or GC labels.
#' @return A ggplot object.
#' @export
read_qc_base_content <- function(data, position = "Position", value = "Ratio",
                                 base = "Base", stage = "Stage",
                                 sample = "SampleID", read = NULL,
                                 colors = NULL, title = "Base content",
                                 style = NULL, output_file = NULL,
                                 width = NULL, height = NULL, dpi = NULL) {
  columns <- c(position, value, base, stage, sample, read)
  columns <- columns[!is.null(columns)]
  .sanshu_require_columns(data, columns)
  .sanshu_numeric(data, c(position, value))
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  if (!is.null(read)) {
    key <- interaction(d[[sample]], d[[stage]], drop = TRUE)
    first <- !grepl("2$|read2|r2", tolower(as.character(d[[read]])))
    maxima <- tapply(d[[position]][first], key[first], max)
    shift <- unname(maxima[as.character(key)])
    second <- !first
    if (any(!is.finite(shift[second]))) {
      stop("Every second-read panel must have matching first-read positions.", call. = FALSE)
    }
    d[[position]][second] <- d[[position]][second] + shift[second]
  }
  bases <- unique(as.character(d[[base]]))
  default <- c(A = "#1F77B4", T = "#FF7F0E", C = "#2CA02C",
               G = "#D62728", N = "#7F7F7F", GC = "#9467BD")
  palette <- .sanshu_palette(bases, colors %||% default, z$style)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[position]], y = .data[[value]], color = .data[[base]],
    group = interaction(.data[[base]], .data[[stage]], .data[[sample]])
  )) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::facet_grid(stats::as.formula(paste(sample, "~", stage)), scales = "free_x") +
    ggplot2::scale_color_manual(values = palette, breaks = bases, drop = FALSE) +
    ggplot2::labs(x = "Position", y = "Base content (%)", color = "Base", title = title) +
    z$style$ggplot_theme
  .sanshu_finish(p, z)
}

#' Draw mapped-read genomic-region composition
#'
#' @param data Data frame containing one row per category.
#' @param category,value Column names for genomic region and percentage/count.
#' @param sample Optional sample column used for faceting.
#' @inheritParams read_qc_error
#' @return A ggplot object.
#' @export
mapping_region_pie <- function(data, category = "Region", value = "Percentage",
                               sample = NULL, colors = NULL,
                               title = "Mapped-read distribution", style = NULL,
                               output_file = NULL, width = NULL, height = NULL,
                               dpi = NULL) {
  .sanshu_require_columns(data, c(category, value, sample)[!is.null(c(category, value, sample))])
  .sanshu_numeric(data, value)
  if (any(data[[value]] < 0)) stop("`value` cannot contain negative values.", call. = FALSE)
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  levels <- unique(as.character(d[[category]]))
  palette <- .sanshu_palette(levels, colors, z$style)
  d$.label <- paste0(d[[category]], " (", sprintf("%.2f", 100 * d[[value]] / ave(
    d[[value]], if (is.null(sample)) rep(1, nrow(d)) else d[[sample]], FUN = sum
  )), "%)")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data[[value]], fill = .data[[category]])) +
    ggplot2::geom_col(width = 1, color = "white", linewidth = 0.3) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$.label),
      position = ggplot2::position_stack(vjust = 0.5),
      family = z$style$global$font_family,
      size = z$style$text$data_label$size / 3.2
    ) +
    ggplot2::scale_fill_manual(values = palette, breaks = levels, guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank())
  if (!is.null(sample)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", sample)))
  .sanshu_finish(p, z)
}

#' Draw gene-body coverage curves
#'
#' @param data Long-form coverage data.
#' @param position,coverage,sample Column names.
#' @inheritParams read_qc_error
#' @return A ggplot object.
#' @export
gene_body_coverage <- function(data, position = "GeneBodyPercent",
                               coverage = "Coverage", sample = "SampleID",
                               colors = NULL, title = "Gene body coverage",
                               style = NULL, output_file = NULL, width = NULL,
                               height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(position, coverage, sample))
  .sanshu_numeric(data, c(position, coverage))
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  levels <- unique(as.character(data[[sample]]))
  palette <- .sanshu_palette(levels, colors, z$style)
  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[[position]], y = .data[[coverage]], color = .data[[sample]],
    group = .data[[sample]]
  )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = palette, breaks = levels, drop = FALSE) +
    ggplot2::labs(x = "Gene body percentile (5' to 3')", y = "Coverage",
                  color = "Sample", title = title) +
    z$style$ggplot_theme
  .sanshu_finish(p, z)
}

#' Draw splice-junction saturation curves
#'
#' @param data Long-form saturation data.
#' @param percent,count,type,sample Column names.
#' @inheritParams read_qc_error
#' @return A ggplot object.
#' @export
junction_saturation <- function(data, percent = "PercentReads",
                                count = "JunctionCount", type = "JunctionType",
                                sample = "SampleID", colors = NULL,
                                title = "Junction saturation", style = NULL,
                                output_file = NULL, width = NULL, height = NULL,
                                dpi = NULL) {
  .sanshu_require_columns(data, c(percent, count, type, sample))
  .sanshu_numeric(data, c(percent, count))
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  levels <- unique(as.character(data[[type]]))
  default <- c("All junctions" = "#0000FF", "known junctions" = "#FF0000",
               "novel junctions" = "#008000")
  palette <- .sanshu_palette(levels, colors %||% default, z$style)
  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[[percent]], y = .data[[count]], color = .data[[type]],
    group = interaction(.data[[type]], .data[[sample]])
  )) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(shape = 1, size = 2) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", sample))) +
    ggplot2::scale_color_manual(values = palette, breaks = levels, drop = FALSE) +
    ggplot2::labs(x = "Percent of total reads",
                  y = "Number of splicing junctions (x1000)",
                  color = NULL, title = title) +
    z$style$ggplot_theme
  .sanshu_finish(p, z)
}

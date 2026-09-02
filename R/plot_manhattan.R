# Alternating chromosome colours for Manhattan plots.
.manhattan_chr_colors <- function(n_chr) {
  rep(c("#3A7CA5", "#6C9FB5"), length.out = n_chr)
}

#' Draw a Manhattan plot
#'
#' Draws a genome-wide association Manhattan plot from a GWAS result table
#' (`CHR`, `BP`, `P`) using [ggplot2]. Chromosome positions are accumulated
#' across chromosomes; two significance thresholds (Bonferroni and suggestive)
#' are annotated by default. A `SNP` column is optional and ignored for the
#' point layer (labels are not drawn to avoid overcrowding).
#'
#' @param data Data frame with the `CHR`, `BP`, and `P` columns.
#' @param chr_column Name of the chromosome column.
#' @param bp_column Name of the base-pair position column.
#' @param p_column Name of the P-value column.
#' @param bonf_threshold P-value for the Bonferroni threshold line (default
#'   `5e-8`). `NULL` hides the line.
#' @param suggest_threshold P-value for the suggestive threshold line (default
#'   `1e-5`). `NULL` hides the line.
#' @param point_alpha Point opacity.
#' @param point_size Point size.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels.
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()].
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object.
#' @export
plot_manhattan <- function(
    data,
    chr_column = "CHR",
    bp_column = "BP",
    p_column = "P",
    bonf_threshold = 5e-8,
    suggest_threshold = 1e-5,
    point_alpha = 0.6,
    point_size = 1.2,
    title = "Manhattan Plot",
    x_label = "Chromosome",
    y_label = NULL,
    font_family = "sans",
    style = NULL,
    output_file = NULL,
    figure_width = 14,
    figure_height = 6,
    dpi = 300) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (column in c(chr_column, bp_column, p_column)) {
    if (!is.character(column) || length(column) != 1L || is.na(column) ||
        !nzchar(column) || !column %in% names(data)) {
      stop("Column `", column, "` is not a column in `data`.", call. = FALSE)
    }
  }
  if (!is.numeric(data[[bp_column]])) {
    stop("`bp_column` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(data[[p_column]])) {
    stop("`p_column` must be numeric.", call. = FALSE)
  }
  .assert_probability(point_alpha, "point_alpha")
  .assert_positive_number(point_size, "point_size")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  if (!is.null(bonf_threshold)) {
    .assert_probability(bonf_threshold, "bonf_threshold")
  }
  if (!is.null(suggest_threshold)) {
    .assert_probability(suggest_threshold, "suggest_threshold")
  }

  data[[chr_column]] <- factor(data[[chr_column]],
                              unique(as.character(data[[chr_column]])))
  chr_levels <- levels(data[[chr_column]])
  chr_colors <- .manhattan_chr_colors(length(chr_levels))

  chr_len <- stats::aggregate(
    data[[bp_column]], by = list(chr = data[[chr_column]]), FUN = max
  )
  chr_tot <- cumsum(as.numeric(chr_len$x)) - as.numeric(chr_len$x)
  names(chr_tot) <- as.character(chr_len$chr)
  data[[".bp_cum"]] <- data[[bp_column]] + chr_tot[as.character(data[[chr_column]])]

  chr_center <- stats::aggregate(
    data[[".bp_cum"]], by = list(chr = data[[chr_column]]),
    FUN = function(x) (min(x) + max(x)) / 2
  )
  chr_breaks <- chr_center$x
  names(chr_breaks) <- as.character(chr_center$chr)

  data[[".minus_log10_p"]] <- -log10(data[[p_column]])

  if (is.null(style)) {
    style <- choose_plot_style(
      font_family = font_family,
      title = list(size = 16),
      axis_title = list(size = 13),
      axis_text = list(size = 11),
      legend = list(show = FALSE),
      panel = list(major_grid = TRUE, minor_grid = FALSE)
    )
  }
  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
  }

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$.bp_cum, y = .data$.minus_log10_p, color = .data[[chr_column]])
  ) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size, shape = 16) +
    ggplot2::scale_color_manual(values = chr_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = unname(chr_breaks), labels = names(chr_breaks),
      expand = c(0.02, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    style$ggplot_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8)
    )

  if (!is.null(bonf_threshold)) {
    y_bonf <- -log10(bonf_threshold)
    plot <- plot +
      ggplot2::geom_hline(yintercept = y_bonf, linetype = "dashed",
                          color = "red", linewidth = 0.7) +
      ggplot2::annotate("text", x = min(data$.bp_cum), y = y_bonf + 0.3,
                        label = sprintf("Bonferroni (%.0e)", bonf_threshold),
                        color = "red", hjust = 0, size = 3)
  }
  if (!is.null(suggest_threshold)) {
    y_suggest <- -log10(suggest_threshold)
    plot <- plot +
      ggplot2::geom_hline(yintercept = y_suggest, linetype = "dashed",
                          color = "blue", linewidth = 0.5) +
      ggplot2::annotate("text", x = min(data$.bp_cum), y = y_suggest + 0.3,
                        label = sprintf("Suggestive (%.0e)", suggest_threshold),
                        color = "blue", hjust = 0, size = 3)
  }

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

#' Draw a QQ plot of observed versus expected P values
#'
#' Draws a Q-Q plot comparing observed to expected -log10 P values and reports
#' the genomic inflation factor (lambda, the median chi-square ratio) in the
#' title.
#'
#' @param data Data frame containing the P-value column.
#' @param p_column Name of the P-value column.
#' @param point_color Point colour.
#' @param point_alpha Point opacity.
#' @param point_size Point size.
#' @param title Optional plot title. If `NULL`, lambda is appended.
#' @param x_label,y_label Axis labels.
#' @param geom_abline Whether to draw the diagonal identity line.
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()].
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object.
#' @export
plot_qq <- function(
    data,
    p_column = "P",
    point_color = "#3A7CA5",
    point_alpha = 0.6,
    point_size = 1.5,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    geom_abline = TRUE,
    font_family = "sans",
    style = NULL,
    output_file = NULL,
    figure_width = 7,
    figure_height = 7,
    dpi = 300) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(p_column) || length(p_column) != 1L || is.na(p_column) ||
      !nzchar(p_column) || !p_column %in% names(data)) {
    stop("Column `", p_column, "` is not a column in `data`.", call. = FALSE)
  }
  if (!is.numeric(data[[p_column]])) {
    stop("`p_column` must be numeric.", call. = FALSE)
  }
  .assert_probability(point_alpha, "point_alpha")
  .assert_positive_number(point_size, "point_size")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  .assert_flag(geom_abline, "geom_abline")
  if (!is.character(point_color) || length(point_color) != 1L ||
      is.na(point_color) || !nzchar(point_color)) {
    stop("`point_color` must be one non-empty colour value.", call. = FALSE)
  }

  p_values <- data[[p_column]]
  p_values <- p_values[is.finite(p_values) & p_values > 0]
  if (!length(p_values)) {
    stop("No positive P values remain for the QQ plot.", call. = FALSE)
  }
  n <- length(p_values)
  expected_p <- (seq_len(n) - 0.5) / n
  qq_data <- data.frame(
    .exp_log10P = -log10(expected_p),
    .obs_log10P = -log10(sort(p_values))
  )

  chisq_obs <- stats::qchisq(1 - p_values, df = 1)
  lambda <- stats::median(chisq_obs) / stats::qchisq(0.5, df = 1)

  if (is.null(title)) {
    title <- sprintf("QQ Plot (lambda = %.4f)", lambda)
  }
  if (is.null(x_label)) {
    x_label <- expression(Expected ~ -log[10](italic(P)))
  }
  if (is.null(y_label)) {
    y_label <- expression(Observed ~ -log[10](italic(P)))
  }

  if (is.null(style)) {
    style <- choose_plot_style(
      font_family = font_family,
      title = list(size = 16),
      axis_title = list(size = 13),
      axis_text = list(size = 11),
      legend = list(show = FALSE)
    )
  }
  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
  }

  plot <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes(x = .data$.exp_log10P, y = .data$.obs_log10P)
  ) +
    ggplot2::geom_point(color = point_color, alpha = point_alpha,
                        size = point_size) +
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    style$ggplot_theme +
    ggplot2::coord_fixed()

  if (geom_abline) {
    plot <- plot + ggplot2::geom_abline(
      intercept = 0, slope = 1, color = "red", linewidth = 0.8,
      linetype = "dashed"
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

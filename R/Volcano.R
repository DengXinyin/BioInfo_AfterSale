#' Plot a differential-expression volcano plot
#'
#' @param data Differential-expression data frame.
#' @param log2fc,pvalue Column names.
#' @param status Optional precomputed status column.
#' @param fc_cutoff,p_cutoff Thresholds used when `status` is NULL.
#' @param status_levels,status_colors Legend order and colors.
#' @param title Plot title; `NULL` hides it.
#' @param legend_title Legend title; `NULL` removes only the title.
#' @param legend_inside Numeric x/y position for the legend inside the panel.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
# All volcano-plot text prefers Times New Roman; the style family controls all labels.
plot_volcano <- function(data, log2fc = "log2FC", pvalue = "pvalue", status = NULL,
                         fc_cutoff = 1, p_cutoff = 0.05,
                         status_levels = c("Up", "Not sig", "Down"),
                         status_colors = c("#E25659", "#CFCFCF", "#335372"),
                         title = NULL, legend_title = NULL, style = NULL,
                         legend_inside = c(0.82, 0.82),
                         output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  style <- .plot_style_or_default(style); d <- as.data.frame(data)
  if (length(miss <- setdiff(c(log2fc, pvalue, status), names(d))))
    stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d$.logp <- -log10(pmax(as.numeric(d[[pvalue]]), .Machine$double.xmin))
  if (is.null(status)) {
    d$.status <- ifelse(d[[log2fc]] >= fc_cutoff & d[[pvalue]] < p_cutoff, "Up",
                        ifelse(d[[log2fc]] <= -fc_cutoff & d[[pvalue]] < p_cutoff, "Down", "Not sig"))
  } else d$.status <- as.character(d[[status]])
  d$.status <- factor(d$.status, levels = status_levels)
  pal <- stats::setNames(status_colors, status_levels)
  counts <- table(factor(d$.status, levels = status_levels))
  legend_labels <- stats::setNames(
    paste0(status_levels, "(", as.integer(counts), ")"), status_levels
  )
  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[log2fc]], .data$.logp, color = .data$.status)) +
    ggplot2::geom_point(size = 1.4, alpha = .75) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggplot2::scale_color_manual(values = pal, labels = legend_labels, drop = FALSE) +
    ggplot2::labs(x = expression(log[2](Fold~Change)), y = expression(-log[10](Pvalue)), title = title, color = legend_title) +
    style$ggplot_theme + ggplot2::theme(
      legend.position = legend_inside,
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank()
    )
  .save_bioinfo_plot(p, style, output_file, width, height, dpi); p
}

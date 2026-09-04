#' Draw a sample-correlation heatmap with group annotations
#'
#' @param correlation Square numeric correlation matrix.
#' @param group Optional named sample-group vector.
#' @param group_colors Optional named group colors.
#' @param cluster Cluster rows and columns together using correlation distance.
#' @param show_values Display correlation coefficients in cells.
#' @param digits Number of displayed decimal places.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
sample_correlation_heatmap <- function(
    correlation, group = NULL, group_colors = NULL, cluster = FALSE,
    show_values = TRUE, digits = 3, title = "Sample correlation",
    style = NULL, output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  matrix <- as.matrix(correlation)
  suppressWarnings(storage.mode(matrix) <- "double")
  if (!nrow(matrix) || nrow(matrix) != ncol(matrix) ||
      is.null(rownames(matrix)) || is.null(colnames(matrix))) {
    stop("`correlation` must be a named square matrix.", call. = FALSE)
  }
  if (!identical(rownames(matrix), colnames(matrix))) {
    stop("Correlation row and column names must be identical and in the same order.", call. = FALSE)
  }
  if (any(!is.finite(matrix)) || any(matrix < -1 | matrix > 1)) {
    stop("Correlation values must be finite and between -1 and 1.", call. = FALSE)
  }
  .assert_flag(cluster, "cluster")
  .assert_flag(show_values, "show_values")
  .assert_positive_number(digits + 1, "digits + 1")
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  samples <- colnames(matrix)
  if (cluster && length(samples) > 2L) {
    order <- stats::hclust(stats::as.dist(1 - matrix))$order
    matrix <- matrix[order, order, drop = FALSE]
    samples <- colnames(matrix)
  }
  if (!is.null(group)) {
    if (is.null(names(group))) stop("`group` must be named by sample.", call. = FALSE)
    missing <- setdiff(samples, names(group))
    extra <- setdiff(names(group), samples)
    if (length(missing) || length(extra)) {
      stop("`group` names must match correlation samples exactly.", call. = FALSE)
    }
    group <- as.character(group[samples])
  }
  d <- as.data.frame(as.table(matrix), stringsAsFactors = FALSE)
  names(d) <- c("Row", "Column", "Correlation")
  d$Row <- factor(d$Row, levels = rev(samples))
  d$Column <- factor(d$Column, levels = samples)
  d$Label <- formatC(d$Correlation, digits = as.integer(digits), format = "f")
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$Column, .data$Row, fill = .data$Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-1, 1), name = "Correlation") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      family = z$style$global$font_family, angle = 90, hjust = 1, vjust = 0.5
    ), panel.grid = ggplot2::element_blank()) +
    ggplot2::coord_equal()
  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$Label),
      family = z$style$global$font_family,
      size = z$style$text$data_label$size / 3.2
    )
  }
  if (!is.null(group)) {
    group_levels <- unique(group)
    palette <- .sanshu_palette(group_levels, group_colors, z$style)
    labels <- paste0(samples, "  [", group, "]")
    p <- p + ggplot2::scale_x_discrete(labels = stats::setNames(labels, samples))
    attr(p, "sample_group_colors") <- palette
  }
  .sanshu_finish(p, z)
}

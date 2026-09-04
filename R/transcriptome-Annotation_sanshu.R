#' Draw a radial annotation-summary flower
#'
#' @param data Annotation database summary.
#' @param database,count Column names for database and annotated-gene count.
#' @param total Core count displayed in the centre.
#' @param colors Optional database colors.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
annotation_flower <- function(data, database = "Database", count = "Count",
                              total, colors = NULL,
                              title = "Functional annotation summary",
                              style = NULL, output_file = NULL, width = NULL,
                              height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(database, count))
  .sanshu_numeric(data, count)
  if (!is.numeric(total) || length(total) != 1L || is.na(total) || total < 0) {
    stop("`total` must be one non-negative number.", call. = FALSE)
  }
  if (any(data[[count]] < 0)) stop("Annotation counts cannot be negative.", call. = FALSE)
  if (anyDuplicated(data[[database]])) stop("Database names must be unique.", call. = FALSE)
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data[order(data[[count]], decreasing = TRUE), , drop = FALSE]
  labels <- as.character(d[[database]])
  palette <- .sanshu_palette(labels, colors, z$style)
  n <- nrow(d)
  theta <- seq(0, 2 * pi, length.out = 121)
  petals <- do.call(rbind, lapply(seq_len(n), function(i) {
    angle <- pi / 2 - 2 * pi * (i - 1) / n
    centre_x <- 2.15 * cos(angle)
    centre_y <- 2.15 * sin(angle)
    major <- 1.35
    minor <- 0.58
    data.frame(
      x = centre_x + major * cos(theta) * cos(angle) - minor * sin(theta) * sin(angle),
      y = centre_y + major * cos(theta) * sin(angle) + minor * sin(theta) * cos(angle),
      Database = labels[i], petal = i
    )
  }))
  text <- data.frame(
    x = 2.15 * cos(pi / 2 - 2 * pi * (seq_len(n) - 1) / n),
    y = 2.15 * sin(pi / 2 - 2 * pi * (seq_len(n) - 1) / n),
    label = format(d[[count]], big.mark = ",", scientific = FALSE)
  )
  outer <- data.frame(
    x = 3.85 * cos(pi / 2 - 2 * pi * (seq_len(n) - 1) / n),
    y = 3.85 * sin(pi / 2 - 2 * pi * (seq_len(n) - 1) / n),
    label = labels
  )
  core <- data.frame(x = 0.95 * cos(theta), y = 0.95 * sin(theta))
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = petals,
      ggplot2::aes(.data$x, .data$y, group = .data$petal, fill = .data$Database),
      alpha = 0.78, color = "white", linewidth = 0.35
    ) +
    ggplot2::geom_polygon(data = core, ggplot2::aes(.data$x, .data$y),
                          fill = "white", color = "#777777", linewidth = 0.6) +
    ggplot2::geom_text(data = text, ggplot2::aes(.data$x, .data$y, label = .data$label),
                       family = z$style$global$font_family,
                       size = z$style$text$data_label$size / 3.2) +
    ggplot2::geom_text(data = outer, ggplot2::aes(.data$x, .data$y, label = .data$label),
                       family = z$style$global$font_family,
                       size = z$style$text$axis_text$size / 3.2) +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = paste0("Core\n", format(total, big.mark = ",", scientific = FALSE)),
      family = z$style$global$font_family,
      size = z$style$text$data_label$size / 3.2
    ) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::coord_equal(xlim = c(-4.7, 4.7), ylim = c(-4.7, 4.7), clip = "off") +
    ggplot2::labs(title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
  .sanshu_finish(p, z)
}

#' Draw a functional-annotation count bar chart
#'
#' @param data Annotation summary table.
#' @param category,value,fill Column names. `fill = NULL` draws one series.
#' @param colors Optional fill colors.
#' @param horizontal Draw horizontal bars.
#' @param position Bar position, either `"stack"` or `"dodge"`.
#' @param top_n Optional number of highest-total categories to display.
#' @inheritParams annotation_flower
#' @return A ggplot object.
#' @export
annotation_bar <- function(data, category = "Category", value = "Count",
                           fill = NULL, colors = NULL, horizontal = TRUE,
                           position = c("stack", "dodge"), top_n = NULL,
                           title = NULL, style = NULL, output_file = NULL,
                           width = NULL, height = NULL, dpi = NULL) {
  position <- match.arg(position)
  .sanshu_require_columns(data, c(category, value, fill))
  .sanshu_numeric(data, value)
  if (any(data[[value]] < 0)) stop("Counts cannot be negative.", call. = FALSE)
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  totals <- tapply(d[[value]], d[[category]], sum)
  ordered <- names(sort(totals, decreasing = TRUE))
  if (!is.null(top_n)) {
    .assert_positive_number(top_n, "top_n")
    ordered <- utils::head(ordered, as.integer(top_n))
    d <- d[as.character(d[[category]]) %in% ordered, , drop = FALSE]
  }
  d[[category]] <- factor(as.character(d[[category]]), levels = rev(ordered))
  if (is.null(fill)) {
    d$.series <- "Count"
    fill <- ".series"
  }
  fill_levels <- unique(as.character(d[[fill]]))
  palette <- .sanshu_palette(fill_levels, colors, z$style)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[category]], y = .data[[value]], fill = .data[[fill]]
  )) +
    ggplot2::geom_col(position = position, width = 0.72) +
    ggplot2::scale_fill_manual(values = palette, breaks = fill_levels,
                               name = if (identical(fill, ".series")) NULL else fill,
                               guide = if (identical(fill, ".series")) "none" else "legend") +
    ggplot2::labs(x = category, y = value, title = title) +
    z$style$ggplot_theme
  if (horizontal) p <- p + ggplot2::coord_flip()
  .sanshu_finish(p, z)
}

#' Draw an annotation composition pie chart
#'
#' @param data Annotation summary table.
#' @param category,value Column names.
#' @param colors Optional category colors.
#' @param top_n Number of largest categories retained; remaining categories are
#'   combined as `other_label`.
#' @param other_label Label for combined lower-frequency categories.
#' @inheritParams annotation_flower
#' @return A ggplot object.
#' @export
annotation_pie <- function(data, category = "Category", value = "Count",
                           colors = NULL, top_n = 10, other_label = "Other",
                           title = NULL, style = NULL, output_file = NULL,
                           width = NULL, height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(category, value))
  .sanshu_numeric(data, value)
  if (any(data[[value]] < 0)) stop("Counts cannot be negative.", call. = FALSE)
  .assert_positive_number(top_n, "top_n")
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  totals <- sort(tapply(data[[value]], data[[category]], sum), decreasing = TRUE)
  keep <- names(utils::head(totals, as.integer(top_n)))
  d <- data.frame(
    Category = c(keep, if (length(totals) > length(keep)) other_label else NULL),
    Count = c(unname(totals[keep]), if (length(totals) > length(keep)) sum(totals[setdiff(names(totals), keep)]) else NULL),
    stringsAsFactors = FALSE
  )
  d$Category <- factor(d$Category, levels = rev(d$Category))
  d$Label <- paste0(as.character(d$Category), " (", sprintf("%.1f", 100 * d$Count / sum(d$Count)), "%)")
  levels <- as.character(d$Category)
  palette <- .sanshu_palette(levels, colors, z$style)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data$Count, fill = .data$Category)) +
    ggplot2::geom_col(width = 1, color = "white", linewidth = 0.35) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = palette, breaks = levels, labels = d$Label,
                               name = NULL, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
  .sanshu_finish(p, z)
}

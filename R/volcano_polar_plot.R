#' Draw an open polar-coordinate volcano plot
#'
#' Maps log2 fold change to angle and adjusted-P-value significance to radius.
#' Each point remains one gene. Down-regulated genes occupy the left sector,
#' non-significant genes remain near the centre, and up-regulated genes occupy
#' the right sector. A gap is retained between the extreme negative and
#' positive directions so opposite effects never meet at a circular seam.
#'
#' @param result Differential-expression data frame.
#' @param gene_column Optional gene-name column retained in `plot$data`.
#' @param log2fc_column,padj_column Fold-change and adjusted-P-value columns.
#' @param log2fc_cutoff,padj_cutoff Differential-expression thresholds.
#' @param fc_limit Absolute fold-change value mapped to the outer angular edge.
#'   Values beyond this limit are visually capped but retained in `plot$data`.
#' @param angular_span Maximum angle in degrees on either side of the centre.
#' @param group_colors Named colors for the three differential groups.
#' @param title Plot title.
#' @param font_family Font family.
#' @param point_size,point_alpha Point appearance.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object containing the classified polar-coordinate data.
#' @export
volcano_polar_plot <- function(
    result,
    gene_column = NULL,
    log2fc_column = "log2FoldChange",
    padj_column = "padj",
    log2fc_cutoff = 1,
    padj_cutoff = 0.05,
    fc_limit = NULL,
    angular_span = 150,
    group_colors = c(
      Up = "#F04438",
      Down = "#3977B8",
      `Not significant` = "#B7B7B7"
    ),
    title = "Polar volcano plot",
    font_family = "sans",
    point_size = 1.25,
    point_alpha = 0.66,
    output_file = NULL,
    figure_width = 8,
    figure_height = 8,
    dpi = 300) {

  if (!is.data.frame(result)) stop("`result` must be a data frame.", call. = FALSE)
  required <- c(log2fc_column, padj_column, gene_column)
  required <- required[!is.na(required) & nzchar(required)]
  missing_columns <- setdiff(required, names(result))
  if (length(missing_columns)) {
    stop("Result is missing column(s): ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  .assert_positive_number(log2fc_cutoff, "log2fc_cutoff")
  .assert_probability(padj_cutoff, "padj_cutoff")
  if (padj_cutoff <= 0) stop("`padj_cutoff` must be greater than 0.", call. = FALSE)
  .assert_positive_number(point_size, "point_size")
  .assert_probability(point_alpha, "point_alpha")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  if (!is.numeric(angular_span) || length(angular_span) != 1L ||
      is.na(angular_span) || angular_span <= 90 || angular_span >= 180) {
    stop("`angular_span` must be between 90 and 180 degrees.", call. = FALSE)
  }

  groups <- c("Up", "Down", "Not significant")
  if (!is.character(group_colors) || is.null(names(group_colors)) ||
      any(!groups %in% names(group_colors))) {
    stop("`group_colors` must contain Up, Down, and Not significant.", call. = FALSE)
  }
  group_colors <- group_colors[groups]

  data <- result
  data$.log2FoldChange <- suppressWarnings(as.numeric(data[[log2fc_column]]))
  data$.padj <- suppressWarnings(as.numeric(data[[padj_column]]))
  valid <- is.finite(data$.log2FoldChange) & is.finite(data$.padj) & data$.padj > 0
  invalid_count <- sum(!valid)
  if (invalid_count) {
    warning(sprintf("%d row(s) with invalid or non-positive padj were excluded.",
                    invalid_count), call. = FALSE)
  }
  data <- data[valid, , drop = FALSE]
  if (!nrow(data)) stop("No valid rows remain.", call. = FALSE)

  data$.minus_log10_padj <- -log10(data$.padj)
  data$.Group <- "Not significant"
  significant <- data$.padj < padj_cutoff
  data$.Group[significant & data$.log2FoldChange > log2fc_cutoff] <- "Up"
  data$.Group[significant & data$.log2FoldChange < -log2fc_cutoff] <- "Down"
  data$.Group <- factor(data$.Group, levels = groups)
  if (!is.null(gene_column)) data$.Gene <- as.character(data[[gene_column]])

  if (is.null(fc_limit)) {
    fc_limit <- max(
      2 * log2fc_cutoff,
      unname(stats::quantile(abs(data$.log2FoldChange), 0.995, na.rm = TRUE))
    )
  }
  .assert_positive_number(fc_limit, "fc_limit")
  if (fc_limit <= log2fc_cutoff) {
    stop("`fc_limit` must exceed `log2fc_cutoff`.", call. = FALSE)
  }

  clipped_fc <- pmax(-fc_limit, pmin(fc_limit, data$.log2FoldChange))
  data$.theta <- clipped_fc / fc_limit * angular_span * pi / 180
  data$.radius <- sqrt(data$.minus_log10_padj)
  data$.polar_x <- data$.radius * sin(data$.theta)
  data$.polar_y <- data$.radius * cos(data$.theta)

  max_significance <- max(data$.minus_log10_padj)
  ring_candidates <- c(-log10(padj_cutoff), 5, 10, 25, 50, 100)
  ring_values <- unique(ring_candidates[ring_candidates <= max_significance])
  if (!length(ring_values)) ring_values <- -log10(padj_cutoff)
  circle_angle <- seq(-pi, pi, length.out = 361L)
  circles <- do.call(rbind, lapply(ring_values, function(value) {
    radius <- sqrt(value)
    data.frame(
      .ring = value,
      x = radius * sin(circle_angle),
      y = radius * cos(circle_angle)
    )
  }))

  ray_fc <- c(-fc_limit, -log2fc_cutoff, 0, log2fc_cutoff, fc_limit)
  ray_theta <- ray_fc / fc_limit * angular_span * pi / 180
  outer_radius <- sqrt(max_significance) * 1.04
  rays <- data.frame(
    x = 0,
    y = 0,
    xend = outer_radius * sin(ray_theta),
    yend = outer_radius * cos(ray_theta),
    label = c(
      sprintf("%.1f", -fc_limit),
      sprintf("-%.1f", log2fc_cutoff),
      "0",
      sprintf("+%.1f", log2fc_cutoff),
      sprintf("+%.1f", fc_limit)
    )
  )
  ray_labels <- transform(
    rays,
    x = outer_radius * 1.08 * sin(ray_theta),
    y = outer_radius * 1.08 * cos(ray_theta)
  )
  ring_labels <- data.frame(
    x = 0,
    y = sqrt(ring_values),
    label = sprintf("%.1f", ring_values)
  )

  # Draw non-significant genes first and highlighted groups on top.
  plot_data <- data[order(data$.Group == "Not significant", decreasing = TRUE), ]
  plot <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = circles,
      ggplot2::aes(x = x, y = y, group = .ring),
      color = "#D9D9D9", linewidth = 0.42
    ) +
    ggplot2::geom_segment(
      data = rays,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "#C9C9C9", linewidth = 0.42
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = .polar_x, y = .polar_y, color = .Group),
      size = point_size, alpha = point_alpha, stroke = 0
    ) +
    ggplot2::geom_text(
      data = ray_labels,
      ggplot2::aes(x = x, y = y, label = label),
      family = font_family, size = 3.2, color = "#444444"
    ) +
    ggplot2::geom_label(
      data = ring_labels,
      ggplot2::aes(x = x, y = y, label = label),
      family = font_family, size = 2.7, color = "#777777",
      fill = "white", linewidth = 0, label.padding = grid::unit(0.08, "lines")
    ) +
    ggplot2::annotate(
      "text", x = 0, y = outer_radius * 1.19,
      label = expression(log[2](FoldChange)~"(angle)"),
      family = font_family, size = 3.5, color = "#333333"
    ) +
    ggplot2::annotate(
      "text", x = 0, y = -outer_radius * 0.22,
      label = expression(-log[10](padj)~"(radius)"),
      family = font_family, size = 3.5, color = "#555555"
    ) +
    ggplot2::scale_color_manual(values = group_colors, breaks = groups, drop = FALSE) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = "Angle = log2FC | Radius = sqrt(-log10 adjusted P value)",
      color = NULL
    ) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10.5, color = "#666666", hjust = 0.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10.5),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(18, 30, 24, 30)
    )

  plot$data <- data
  attr(plot, "fc_limit") <- fc_limit
  attr(plot, "excluded_rows") <- invalid_count
  if (!is.null(output_file)) {
    .volcano_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

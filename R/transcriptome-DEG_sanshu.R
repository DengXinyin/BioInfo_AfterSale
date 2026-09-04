#' Draw differential-expression counts by comparison
#'
#' @param data Count table.
#' @param comparison,direction,count Column names.
#' @param colors Optional direction colors.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
deg_count_bar <- function(data, comparison = "Comparison", direction = "Direction",
                          count = "Count", colors = c(Up = "#D5695D", Down = "#65A479"),
                          title = "Differentially expressed genes", style = NULL,
                          output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(comparison, direction, count))
  .sanshu_numeric(data, count)
  if (any(data[[count]] < 0)) stop("DEG counts cannot be negative.", call. = FALSE)
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  totals <- sort(tapply(d[[count]], d[[comparison]], sum), decreasing = TRUE)
  d[[comparison]] <- factor(as.character(d[[comparison]]), levels = names(totals))
  directions <- unique(as.character(d[[direction]]))
  palette <- .sanshu_palette(directions, colors, z$style)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[comparison]], y = .data[[count]], fill = .data[[direction]]
  )) +
    ggplot2::geom_col(position = "dodge", width = 0.72) +
    ggplot2::scale_fill_manual(values = palette, breaks = directions, name = NULL) +
    ggplot2::labs(x = "Comparison", y = "Gene count", title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      family = z$style$global$font_family, angle = 45, hjust = 1
    ))
  .sanshu_finish(p, z)
}

#' Draw a DEG comparison flower for many sets
#'
#' @param sets Named list of DEG identifiers.
#' @param colors Optional set colors.
#' @inheritParams deg_count_bar
#' @return A ggplot object.
#' @export
deg_flower <- function(sets, colors = NULL, title = "DEG set flower",
                       style = NULL, output_file = NULL, width = NULL,
                       height = NULL, dpi = NULL) {
  if (!is.list(sets) || length(sets) < 3L || is.null(names(sets)) ||
      any(!nzchar(names(sets)))) {
    stop("`sets` must be a named list containing at least three DEG sets.", call. = FALSE)
  }
  sets <- lapply(sets, function(x) unique(as.character(x[!is.na(x)])))
  summary <- data.frame(Comparison = names(sets), Count = lengths(sets))
  core <- length(Reduce(intersect, sets))
  annotation_flower(
    summary, database = "Comparison", count = "Count", total = core,
    colors = colors, title = title, style = style, output_file = output_file,
    width = width, height = height, dpi = dpi
  )
}

#' Draw a differential-expression MA plot
#'
#' @param data Differential-expression table.
#' @param mean_expression,log2fc,status Column names.
#' @param colors Optional status colors.
#' @param fc_cutoff Fold-change guide-line threshold.
#' @inheritParams deg_count_bar
#' @return A ggplot object.
#' @export
deg_ma <- function(data, mean_expression = "MeanExpression", log2fc = "log2FC",
                   status = "Status",
                   colors = c(Up = "#D5695D", NoSig = "#BDBDBD", Down = "#65A479"),
                   fc_cutoff = 1, title = "MA plot", style = NULL,
                   output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(mean_expression, log2fc, status))
  .sanshu_numeric(data, c(mean_expression, log2fc))
  if (any(data[[mean_expression]] <= 0)) {
    stop("Mean expression values must be positive before log2 transformation.", call. = FALSE)
  }
  .assert_positive_number(fc_cutoff, "fc_cutoff")
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  d$.log_mean <- log2(d[[mean_expression]])
  levels <- unique(as.character(d[[status]]))
  palette <- .sanshu_palette(levels, colors, z$style)
  counts <- table(d[[status]])
  labels <- paste0(names(counts), ": ", as.integer(counts))
  names(labels) <- names(counts)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data$.log_mean, y = .data[[log2fc]], color = .data[[status]]
  )) +
    ggplot2::geom_point(size = 0.8, alpha = 0.55) +
    ggplot2::geom_hline(yintercept = c(-fc_cutoff, fc_cutoff),
                        color = "#777777", linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = palette, breaks = levels,
                                labels = labels[levels], name = NULL) +
    ggplot2::labs(x = "log2(Mean of normalized counts)", y = "log2(FoldChange)",
                  title = title) +
    z$style$ggplot_theme
  .sanshu_finish(p, z)
}

#' Draw an UpSet-style intersection summary without extra dependencies
#'
#' @param sets Named list of feature vectors.
#' @param max_intersections Maximum number of non-empty intersections shown.
#' @inheritParams deg_count_bar
#' @return A ggplot object.
#' @export
deg_upset <- function(sets, max_intersections = 30, title = "DEG intersections",
                      style = NULL, output_file = NULL, width = NULL,
                      height = NULL, dpi = NULL) {
  if (!is.list(sets) || length(sets) < 2L || is.null(names(sets)) ||
      any(!nzchar(names(sets)))) {
    stop("`sets` must be a named list containing at least two sets.", call. = FALSE)
  }
  .assert_positive_number(max_intersections, "max_intersections")
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  sets <- lapply(sets, function(x) unique(as.character(x[!is.na(x)])))
  universe <- unique(unlist(sets, use.names = FALSE))
  membership <- vapply(sets, function(x) universe %in% x, logical(length(universe)))
  keys <- apply(membership, 1, function(x) paste0(as.integer(x), collapse = ""))
  counts <- sort(table(keys), decreasing = TRUE)
  counts <- counts[names(counts) != paste0(rep("0", length(sets)), collapse = "")]
  counts <- utils::head(counts, as.integer(max_intersections))
  bits <- do.call(rbind, strsplit(names(counts), "", fixed = TRUE)) == "1"
  n_set <- length(sets)
  n_int <- length(counts)
  matrix_data <- expand.grid(Intersection = seq_len(n_int), Set = seq_len(n_set))
  matrix_data$Active <- as.vector(bits)
  active <- matrix_data[matrix_data$Active, , drop = FALSE]
  spans <- do.call(rbind, lapply(seq_len(n_int), function(i) {
    y <- which(bits[i, ])
    if (length(y) < 2L) return(NULL)
    data.frame(Intersection = i, ymin = min(y), ymax = max(y))
  }))
  maximum <- max(as.numeric(counts))
  bars <- data.frame(
    Intersection = seq_len(n_int), ymin = n_set + 1,
    ymax = n_set + 1 + 4 * as.numeric(counts) / maximum,
    Count = as.numeric(counts)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = matrix_data, ggplot2::aes(.data$Intersection, .data$Set),
                        color = "#D9D9D9", size = 2.6) +
    ggplot2::geom_segment(data = spans, ggplot2::aes(
      x = .data$Intersection, xend = .data$Intersection,
      y = .data$ymin, yend = .data$ymax
    ), linewidth = 0.8) +
    ggplot2::geom_point(data = active, ggplot2::aes(.data$Intersection, .data$Set),
                        color = "black", size = 2.6) +
    ggplot2::geom_rect(data = bars, ggplot2::aes(
      xmin = .data$Intersection - 0.35, xmax = .data$Intersection + 0.35,
      ymin = .data$ymin, ymax = .data$ymax
    ), fill = "black") +
    ggplot2::geom_text(data = bars, ggplot2::aes(
      x = .data$Intersection, y = .data$ymax, label = .data$Count
    ), vjust = -0.25, family = z$style$global$font_family,
    size = z$style$text$data_label$size / 3.2) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_set), labels = names(sets),
      limits = c(0.5, n_set + 5.7), expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_int)) +
    ggplot2::labs(x = "Intersection", y = NULL, title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  .sanshu_finish(p, z)
}

#' Draw clustered gene-expression trends
#'
#' @param data Long-form expression table.
#' @param gene,sample,expression,cluster Column names.
#' @param sample_order Explicit sample order; must match all sample labels.
#' @param line_color,center_color Colors for gene and cluster-centre lines.
#' @inheritParams deg_count_bar
#' @return A ggplot object.
#' @export
expression_trend <- function(data, gene = "Gene", sample = "Sample",
                             expression = "Zscore", cluster = "Cluster",
                             sample_order = NULL, line_color = "#BDBDBD",
                             center_color = "#D5695D", title = "Expression trends",
                             style = NULL, output_file = NULL, width = NULL,
                             height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(gene, sample, expression, cluster))
  .sanshu_numeric(data, expression)
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  observed <- unique(as.character(d[[sample]]))
  if (is.null(sample_order)) sample_order <- observed
  if (!setequal(sample_order, observed) || anyDuplicated(sample_order)) {
    stop("`sample_order` must contain every sample exactly once.", call. = FALSE)
  }
  d[[sample]] <- factor(as.character(d[[sample]]), levels = sample_order)
  centre <- stats::aggregate(d[[expression]], d[c(cluster, sample)], mean)
  names(centre)[ncol(centre)] <- ".centre"
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[sample]], y = .data[[expression]], group = .data[[gene]]
  )) +
    ggplot2::geom_line(color = line_color, alpha = 0.35, linewidth = 0.35) +
    ggplot2::geom_line(data = centre, ggplot2::aes(
      x = .data[[sample]], y = .data$.centre, group = 1
    ), inherit.aes = FALSE, color = center_color, linewidth = 1.1) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", cluster)), scales = "free_y") +
    ggplot2::labs(x = "Sample", y = "Gene expression (Z-score)", title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      family = z$style$global$font_family, angle = 45, hjust = 1
    ))
  .sanshu_finish(p, z)
}

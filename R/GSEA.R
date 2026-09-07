#' Calculate a GSEA running enrichment score
#'
#' Reconstructs the running-score curve from a named ranked metric and one
#' gene set. This is useful when an archived workflow provides an `fgsea`
#' table or plain vectors rather than a `clusterProfiler::gseaResult` object.
#'
#' @param ranked_metric Named numeric vector. Values are sorted decreasingly
#'   unless `sort = FALSE`.
#' @param gene_set Character vector containing members of one gene set.
#' @param exponent Non-negative weight exponent. `1` is weighted GSEA and `0`
#'   gives the unweighted statistic.
#' @param sort Sort `ranked_metric` decreasingly before calculation.
#' @param na_rm Remove non-finite ranked values.
#'
#' @return A `bioinfo_gsea_curve` list containing the ordered metric, matched
#'   gene set, hit positions, running score, peak position/ES, and leading edge.
#' @export
calculate_gsea_curve <- function(ranked_metric, gene_set, exponent = 1,
                                 sort = TRUE, na_rm = TRUE) {
  if (!is.numeric(ranked_metric) || !length(ranked_metric) ||
      is.null(names(ranked_metric)) || any(!nzchar(names(ranked_metric)))) {
    stop("`ranked_metric` must be a non-empty named numeric vector.", call. = FALSE)
  }
  if (anyDuplicated(names(ranked_metric))) {
    stop("`ranked_metric` gene names must be unique.", call. = FALSE)
  }
  if (length(exponent) != 1L || !is.numeric(exponent) || !is.finite(exponent) ||
      exponent < 0) {
    stop("`exponent` must be one non-negative number.", call. = FALSE)
  }
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort) ||
      !is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`sort` and `na_rm` must each be TRUE or FALSE.", call. = FALSE)
  }
  finite <- is.finite(ranked_metric)
  if (!all(finite)) {
    if (!na_rm) stop("`ranked_metric` contains non-finite values.", call. = FALSE)
    ranked_metric <- ranked_metric[finite]
  }
  if (sort) ranked_metric <- base::sort(ranked_metric, decreasing = TRUE)
  gene_set <- unique(stats::na.omit(as.character(gene_set)))
  gene_set <- intersect(gene_set[nzchar(gene_set)], names(ranked_metric))
  n <- length(ranked_metric)
  hits <- which(names(ranked_metric) %in% gene_set)
  if (!length(hits)) stop("`gene_set` has no member in `ranked_metric`.", call. = FALSE)
  if (length(hits) == n) stop("`gene_set` cannot contain every ranked gene.", call. = FALSE)

  in_set <- seq_len(n) %in% hits
  hit_weights <- abs(ranked_metric[in_set])^exponent
  if (!is.finite(sum(hit_weights)) || sum(hit_weights) <= 0) {
    hit_weights <- rep(1, length(hits))
  }
  steps <- rep(-1 / (n - length(hits)), n)
  steps[in_set] <- hit_weights / sum(hit_weights)
  running_score <- cumsum(steps)
  running_score[n] <- 0
  peak_position <- which.max(abs(running_score))
  enrichment_score <- running_score[peak_position]
  leading_hits <- if (enrichment_score >= 0) hits[hits <= peak_position] else hits[hits >= peak_position]

  structure(
    list(
      ranked_metric = ranked_metric,
      gene_set = names(ranked_metric)[hits],
      hits = hits,
      running_score = running_score,
      peak_position = peak_position,
      enrichment_score = enrichment_score,
      leading_edge = names(ranked_metric)[leading_hits],
      exponent = exponent
    ),
    class = c("bioinfo_gsea_curve", "list")
  )
}

#' Plot a three-panel GSEA running-score figure
#'
#' @param running_score Numeric running enrichment score.
#' @param hits Integer hit positions or a logical vector.
#' @param ranked_metric Numeric ranked metric of the same length.
#' @param title Optional title.
#' @param statistics Optional named values displayed only in the top panel.
#' @param colors Positive and negative colors.
#' @param metric_colors Low, midpoint, and high colors for the ranked-metric
#'   band in the middle panel.
#' @param show_metric_band Draw the ranked-metric color band behind hit marks.
#' @param peak_position Optional ES peak rank marked by a dashed line.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional PDF or PNG output filename.
#' @param width,height,dpi Optional output overrides.
#'
#' @return A faceted ggplot object.
#' @export
plot_gsea <- function(running_score, hits, ranked_metric, title = NULL,
                      statistics = NULL, colors = c("#E25659", "#335372"),
                      metric_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                      show_metric_band = TRUE, peak_position = NULL,
                      style = NULL, output_file = NULL, width = NULL,
                      height = NULL, dpi = NULL) {
  style <- .plot_style_or_default(style)
  running_score <- as.numeric(running_score)
  ranked_metric <- as.numeric(ranked_metric)
  n <- length(running_score)
  if (!n || length(ranked_metric) != n) {
    stop("`running_score` and `ranked_metric` must have equal positive lengths.", call. = FALSE)
  }
  if (any(!is.finite(running_score)) || any(!is.finite(ranked_metric))) {
    stop("GSEA vectors must contain only finite values.", call. = FALSE)
  }
  if (is.logical(hits)) {
    if (length(hits) != n) stop("Logical `hits` must match the vector length.", call. = FALSE)
    hits <- which(hits)
  }
  hits <- unique(as.integer(hits))
  if (!length(hits) || anyNA(hits) || any(hits < 1 | hits > n)) {
    stop("`hits` must contain valid ranked positions.", call. = FALSE)
  }
  if (!is.null(peak_position) &&
      (length(peak_position) != 1L || !is.finite(peak_position) ||
       peak_position < 1 || peak_position > n)) {
    stop("`peak_position` must be one valid ranked position.", call. = FALSE)
  }
  if (length(colors) != 2L || length(metric_colors) != 3L) {
    stop("`colors` needs two colors and `metric_colors` needs three.", call. = FALSE)
  }
  panels <- factor(
    c(rep("Running Enrichment Score", n), rep("Hits", length(hits)),
      rep("Ranked Metric", n)),
    levels = c("Running Enrichment Score", "Hits", "Ranked Metric")
  )
  d <- data.frame(
    x = c(seq_len(n), hits, seq_len(n)),
    y = c(running_score, rep(0.86, length(hits)), ranked_metric),
    panel = panels
  )
  band <- data.frame(
    x = seq_len(n), y = rep(0.32, n), metric = ranked_metric,
    panel = factor("Hits", levels = levels(panels))
  )

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y))
  if (show_metric_band) {
    max_abs <- max(abs(ranked_metric))
    if (max_abs == 0) max_abs <- 1
    p <- p + ggplot2::geom_tile(
      data = band,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$metric),
      inherit.aes = FALSE, width = 1.05, height = 0.58
    ) + ggplot2::scale_fill_gradientn(
      colors = metric_colors, values = c(0, 0.5, 1),
      limits = c(-max_abs, max_abs), guide = "none"
    )
  }
  p <- p +
    ggplot2::geom_line(
      data = d[d$panel == "Running Enrichment Score", ],
      color = colors[1], linewidth = 0.9
    ) +
    ggplot2::geom_hline(
      data = data.frame(
        yintercept = 0,
        panel = factor(c("Running Enrichment Score", "Ranked Metric"), levels = levels(panels))
      ),
      ggplot2::aes(yintercept = .data$yintercept),
      colour = "grey55", linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = d[d$panel == "Hits", ],
      ggplot2::aes(xend = .data$x, y = 0.68, yend = 1),
      color = "black", linewidth = 0.28
    ) +
    ggplot2::geom_area(
      data = d[d$panel == "Ranked Metric" & d$y >= 0, ],
      fill = colors[1], alpha = 0.85
    ) +
    ggplot2::geom_area(
      data = d[d$panel == "Ranked Metric" & d$y < 0, ],
      fill = colors[2], alpha = 0.85
    )
  if (!is.null(peak_position)) {
    peak_data <- data.frame(
      xintercept = peak_position,
      panel = factor(c("Running Enrichment Score", "Hits"), levels = levels(panels))
    )
    p <- p + ggplot2::geom_vline(
      data = peak_data, ggplot2::aes(xintercept = .data$xintercept),
      colour = colors[1], linetype = "dashed", linewidth = 0.45
    )
  }
  p <- p +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$panel), scales = "free_y", switch = "y",
      labeller = ggplot2::labeller(panel = c(
        `Running Enrichment Score` = "Running ES",
        Hits = "Hits",
        `Ranked Metric` = "Ranked metric"
      ))
    ) +
    ggplot2::labs(x = "Rank in Ordered Gene List", y = NULL, title = title) +
    style$ggplot_theme +
    ggplot2::theme(
      panel.spacing.y = grid::unit(0, "pt"),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(
        family = style$global$font_family,
        size = min(style$text$axis_text$size, 10),
        face = "plain"
      )
    )
  if (!is.null(statistics)) {
    if (is.null(names(statistics)) || any(!nzchar(names(statistics)))) {
      stop("`statistics` must be a named vector.", call. = FALSE)
    }
    label <- paste(names(statistics), format(statistics, digits = 3), sep = " = ", collapse = "\n")
    stat_data <- data.frame(
      x = Inf, y = Inf, label = label,
      panel = factor("Running Enrichment Score", levels = levels(panels))
    )
    p <- p + ggplot2::geom_text(
      data = stat_data,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      inherit.aes = FALSE, hjust = 1.05, vjust = 1.15,
      family = style$global$font_family,
      size = style$text$data_label$size / 3.2
    )
  }
  attr(p, "peak_position") <- peak_position
  .save_bioinfo_plot(p, style, output_file, width, height, dpi)
  p
}

#' Plot one pathway from ranked values and a gene set
#'
#' @inheritParams calculate_gsea_curve
#' @param title,statistics,colors,metric_colors,show_metric_band,style,output_file,width,height,dpi
#'   Passed to [plot_gsea()].
#'
#' @return A ggplot object with its calculated curve in the `gsea_curve`
#'   attribute.
#' @export
plot_gsea_pathway <- function(ranked_metric, gene_set, exponent = 1,
                              title = NULL, statistics = NULL,
                              colors = c("#E25659", "#335372"),
                              metric_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                              show_metric_band = TRUE, style = NULL,
                              output_file = NULL, width = NULL, height = NULL,
                              dpi = NULL) {
  curve <- calculate_gsea_curve(ranked_metric, gene_set, exponent = exponent)
  p <- plot_gsea(
    running_score = curve$running_score,
    hits = curve$hits,
    ranked_metric = curve$ranked_metric,
    title = title,
    statistics = statistics,
    colors = colors,
    metric_colors = metric_colors,
    show_metric_band = show_metric_band,
    peak_position = curve$peak_position,
    style = style,
    output_file = output_file,
    width = width,
    height = height,
    dpi = dpi
  )
  attr(p, "gsea_curve") <- curve
  p
}

.gsea_result_components <- function(result, gene_set_id) {
  if (!methods::is(result, "gseaResult")) {
    stop("`result` must be a clusterProfiler `gseaResult` object.", call. = FALSE)
  }
  table <- methods::slot(result, "result")
  gene_list <- methods::slot(result, "geneList")
  gene_sets <- methods::slot(result, "geneSets")
  id_match <- which(as.character(table$ID) == gene_set_id)
  if (!length(id_match) && "Description" %in% names(table)) {
    id_match <- which(as.character(table$Description) == gene_set_id)
  }
  if (length(id_match) != 1L) {
    stop("`gene_set_id` must uniquely match an ID or Description in `result`.", call. = FALSE)
  }
  row <- table[id_match, , drop = FALSE]
  id <- as.character(row$ID[1])
  if (!id %in% names(gene_sets)) stop("The selected gene set is absent from `result@geneSets`.", call. = FALSE)
  stats <- c(
    NES = as.numeric(row$NES[1]),
    `P value` = as.numeric(row$pvalue[1]),
    `Adjusted P` = as.numeric(row$p.adjust[1])
  )
  list(id = id, title = as.character(row$Description[1]), statistics = stats,
       ranked_metric = gene_list, gene_set = gene_sets[[id]])
}

#' Plot a pathway stored in a clusterProfiler GSEA result
#'
#' @param result A `clusterProfiler::gseaResult` object.
#' @param gene_set_id One pathway ID or exact Description.
#' @param title Optional title replacing the result Description.
#' @param ... Additional arguments passed to [plot_gsea_pathway()].
#'
#' @return A ggplot object.
#' @export
plot_gsea_result <- function(result, gene_set_id, title = NULL, ...) {
  z <- .gsea_result_components(result, gene_set_id)
  plot_gsea_pathway(
    ranked_metric = z$ranked_metric,
    gene_set = z$gene_set,
    title = title %||% z$title,
    statistics = z$statistics,
    ...
  )
}

#' Batch-export selected pathways from a GSEA result
#'
#' @param result A `clusterProfiler::gseaResult` object.
#' @param gene_set_ids Character vector of pathway IDs or exact Descriptions.
#' @param output_dir Existing or new output directory.
#' @param format Output format: `"pdf"` or `"png"`.
#' @param filename_prefix Filename prefix.
#' @param width,height,dpi Output dimensions and PNG resolution.
#' @param ... Additional arguments passed to [plot_gsea_result()].
#'
#' @return Invisibly, a named character vector of output files.
#' @export
plot_gsea_batch <- function(result, gene_set_ids, output_dir,
                            format = c("pdf", "png"), filename_prefix = "GSEA",
                            width = 9.2, height = 6.8, dpi = 300, ...) {
  format <- match.arg(format)
  if (!length(gene_set_ids) || anyNA(gene_set_ids) || any(!nzchar(gene_set_ids))) {
    stop("`gene_set_ids` must contain at least one non-empty ID.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- vapply(gene_set_ids, function(id) {
    z <- .gsea_result_components(result, id)
    safe <- gsub("[^A-Za-z0-9_.-]+", "_", paste(filename_prefix, z$id, z$title, sep = "_"))
    file <- file.path(output_dir, paste0(safe, ".", format))
    plot_gsea_result(result, id, output_file = file, width = width,
                     height = height, dpi = dpi, ...)
    normalizePath(file, mustWork = TRUE)
  }, character(1))
  names(files) <- gene_set_ids
  invisible(files)
}

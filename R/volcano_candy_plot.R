#' Draw a candy-shaped, group-guided UMAP volcano variant
#'
#' Builds a two-dimensional UMAP embedding from differential-expression
#' features and arranges the three volcano classes as a candy-like composition:
#' down-regulated genes on the left, non-significant genes in the centre, and
#' up-regulated genes on the right. This is an exploratory, group-guided
#' visualization and does not replace a conventional volcano plot.
#'
#' Rows with missing, non-finite, or non-positive adjusted P values are removed
#' before transformation. Gene names are retained as identifiers and labels;
#' they are not numerically encoded as a distance feature.
#'
#' @param result Data frame containing differential-expression results.
#' @param gene_column Optional gene-name column used for labels.
#' @param log2fc_column,padj_column Fold-change and adjusted-P-value columns.
#' @param base_mean_column Optional mean-expression column added as a feature
#'   after a `log1p` transformation.
#' @param log2fc_cutoff,padj_cutoff Differential-expression thresholds.
#' @param n_neighbors UMAP neighbourhood size.
#' @param min_dist UMAP minimum distance.
#' @param target_weight Strength of supervised group guidance in UMAP.
#' @param candy_strength Strength of the final left-centre-right arrangement.
#'   Zero keeps raw UMAP coordinates; one gives the clearest candy silhouette.
#' @param umap_columns Optional pair of precomputed UMAP-coordinate column
#'   names. When supplied, UMAP is not recalculated in R. This supports a
#'   Python-compute/R-plot workflow.
#' @param seed Random seed for reproducible UMAP coordinates.
#' @param label_top Number of most significant genes to label per Up/Down group.
#' @param group_colors Named colors for `Up`, `Down`, and `Not significant`.
#' @param title Plot title.
#' @param font_family Font family used by the plot theme.
#' @param point_size,point_alpha Point appearance.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object. The filtered input and embedding are in `plot$data`.
#' @export
volcano_candy_plot <- function(
    result,
    gene_column = NULL,
    log2fc_column = "log2FoldChange",
    padj_column = "padj",
    base_mean_column = NULL,
    log2fc_cutoff = 1,
    padj_cutoff = 0.05,
    n_neighbors = 35,
    min_dist = 0.18,
    target_weight = 0.65,
    candy_strength = 0.82,
    umap_columns = NULL,
    seed = 20260825,
    label_top = 5,
    group_colors = c(
      Up = "#F04438",
      Down = "#3977B8",
      `Not significant` = "#A9A9A9"
    ),
    title = "Candy volcano: group-guided UMAP",
    font_family = "sans",
    point_size = 1.7,
    point_alpha = 0.72,
    output_file = NULL,
    figure_width = 9,
    figure_height = 6.5,
    dpi = 300) {

  if (!is.data.frame(result)) stop("`result` must be a data frame.", call. = FALSE)
  if (is.null(umap_columns) && !requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for volcano_candy_plot().", call. = FALSE)
  }
  if (!is.null(umap_columns) &&
      (!is.character(umap_columns) || length(umap_columns) != 2L ||
       anyNA(umap_columns) || any(!nzchar(umap_columns)))) {
    stop("`umap_columns` must be NULL or two non-empty column names.", call. = FALSE)
  }
  required <- c(
    log2fc_column, padj_column, gene_column, base_mean_column, umap_columns
  )
  required <- required[!is.na(required) & nzchar(required)]
  missing_columns <- setdiff(required, names(result))
  if (length(missing_columns)) {
    stop("Result is missing column(s): ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  .assert_positive_number(log2fc_cutoff, "log2fc_cutoff")
  .assert_probability(padj_cutoff, "padj_cutoff")
  if (padj_cutoff <= 0) stop("`padj_cutoff` must be greater than 0.", call. = FALSE)
  .assert_probability(min_dist, "min_dist")
  .assert_probability(target_weight, "target_weight")
  .assert_probability(candy_strength, "candy_strength")
  .assert_positive_number(n_neighbors, "n_neighbors")
  .assert_positive_number(point_size, "point_size")
  .assert_probability(point_alpha, "point_alpha")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  if (length(label_top) != 1L || is.na(label_top) || label_top < 0) {
    stop("`label_top` must be one non-negative number.", call. = FALSE)
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
  if (!is.null(base_mean_column)) {
    data$.baseMean <- suppressWarnings(as.numeric(data[[base_mean_column]]))
    valid <- valid & is.finite(data$.baseMean) & data$.baseMean >= 0
  }
  invalid_count <- sum(!valid)
  if (invalid_count) {
    warning(sprintf("%d row(s) with invalid or non-positive values were excluded.",
                    invalid_count), call. = FALSE)
  }
  data <- data[valid, , drop = FALSE]
  if (nrow(data) < 10L) stop("At least 10 valid rows are required.", call. = FALSE)

  data$.minus_log10_padj <- -log10(data$.padj)
  data$.Group <- "Not significant"
  significant <- data$.padj < padj_cutoff
  data$.Group[significant & data$.log2FoldChange >= log2fc_cutoff] <- "Up"
  data$.Group[significant & data$.log2FoldChange <= -log2fc_cutoff] <- "Down"
  data$.Group <- factor(data$.Group, levels = groups)
  if (any(table(data$.Group) == 0L)) {
    stop("All three groups need at least one gene for the candy layout.", call. = FALSE)
  }

  feature_data <- data.frame(
    log2_fc = data$.log2FoldChange,
    abs_log2_fc = abs(data$.log2FoldChange),
    significance = data$.minus_log10_padj,
    signed_significance = sign(data$.log2FoldChange) * data$.minus_log10_padj
  )
  if (!is.null(base_mean_column)) feature_data$log_expression <- log1p(data$.baseMean)
  feature_matrix <- scale(as.matrix(feature_data))

  if (is.null(umap_columns)) {
    neighbours <- max(2L, min(as.integer(n_neighbors), nrow(data) - 1L))
    embedding <- uwot::umap(
      feature_matrix,
      n_neighbors = neighbours,
      min_dist = min_dist,
      metric = "euclidean",
      y = data$.Group,
      target_metric = "categorical",
      target_weight = target_weight,
      n_threads = 1,
      seed = seed,
      verbose = FALSE
    )
    embedding <- scale(embedding)

    # Preserve local UMAP structure while pulling group centroids into a stable
    # left-centre-right candy arrangement.
    anchors <- c(Up = 2.7, Down = -2.7, `Not significant` = 0)
    raw_x <- embedding[, 1]
    raw_y <- embedding[, 2]
    data$.UMAP_raw_1 <- raw_x
    data$.UMAP_raw_2 <- raw_y
    centred_x <- ave(raw_x, data$.Group, FUN = function(x) x - mean(x))
    centred_y <- ave(raw_y, data$.Group, FUN = function(x) x - mean(x))
    candy_x <- anchors[as.character(data$.Group)] + centred_x * 0.72
    candy_y <- centred_y * ifelse(data$.Group == "Not significant", 0.72, 1.0)
    data$.UMAP_1 <- (1 - candy_strength) * raw_x + candy_strength * candy_x
    data$.UMAP_2 <- (1 - candy_strength) * raw_y + candy_strength * candy_y
  } else {
    data$.UMAP_1 <- suppressWarnings(as.numeric(data[[umap_columns[[1]]]]))
    data$.UMAP_2 <- suppressWarnings(as.numeric(data[[umap_columns[[2]]]]))
    if (any(!is.finite(data$.UMAP_1)) || any(!is.finite(data$.UMAP_2))) {
      stop("Precomputed UMAP columns must contain only finite numbers.", call. = FALSE)
    }
    data$.UMAP_raw_1 <- data$.UMAP_1
    data$.UMAP_raw_2 <- data$.UMAP_2
  }

  if (is.null(gene_column)) {
    data$.GeneLabel <- rep("", nrow(data))
  } else {
    data$.GeneLabel <- as.character(data[[gene_column]])
  }
  data$.ShowLabel <- ""
  if (!is.null(gene_column) && label_top > 0) {
    label_rows <- unlist(lapply(c("Up", "Down"), function(group) {
      candidates <- which(data$.Group == group)
      candidates[head(order(data$.padj[candidates]), as.integer(label_top))]
    }), use.names = FALSE)
    data$.ShowLabel[label_rows] <- data$.GeneLabel[label_rows]
  }

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .UMAP_1, y = .UMAP_2, color = .Group, fill = .Group)
  ) +
    ggplot2::stat_ellipse(
      geom = "polygon", type = "norm", level = 0.82,
      alpha = 0.09, linewidth = 0.55, show.legend = FALSE
    ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
    ggplot2::scale_color_manual(values = group_colors, breaks = groups, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_colors, breaks = groups, drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = "Group-guided UMAP of differential-expression features",
      x = "UMAP 1",
      y = "UMAP 2",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, color = "#666666", hjust = 0.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(18, 24, 18, 24)
    )

  if (!is.null(gene_column) && label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    plot <- plot + ggrepel::geom_text_repel(
      ggplot2::aes(label = .ShowLabel),
      size = 3.1, family = font_family, color = "#333333",
      min.segment.length = 0, box.padding = 0.35, point.padding = 0.2,
      max.overlaps = Inf, seed = seed, show.legend = FALSE
    )
  }

  attr(plot, "embedding_method") <- if (is.null(umap_columns)) {
    "group-guided UMAP calculated in R"
  } else {
    "precomputed UMAP coordinates"
  }
  attr(plot, "feature_columns") <- names(feature_data)
  attr(plot, "excluded_rows") <- invalid_count
  if (!is.null(output_file)) {
    .volcano_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

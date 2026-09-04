# Toolbox group palette (matches 组学可视化百宝箱 colour conventions).
.toolbox_group_palette <- function() {
  c(
    "#4DBBD5", "#E64B35", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
  )
}

#' Draw a principal-component-analysis score plot
#'
#' Performs PCA with [stats::prcomp()] on a user-supplied expression matrix
#' and draws a two-component score plot with optional per-group 95% confidence
#' ellipses and optional sample labels. The colour scheme follows the shared
#' toolbox palette used across the 组学可视化百宝箱 modules by default.
#'
#' Columns identified as the sample and group columns are held out of the
#' PCA matrix. When `expr_columns` is `NULL`, every remaining numeric column is
#' used. Rows with any non-finite expression value are removed before scaling.
#'
#' @param data Data frame containing the expression matrix plus, optionally,
#'   sample and group columns.
#' @param sample_column Name of the sample identifier column.
#' @param group_column Name of the grouping column used to colour points.
#' @param expr_columns Optional character vector of expression columns to feed
#'   into [stats::prcomp()]. `NULL` uses every numeric column other than the
#'   sample and group columns.
#' @param center,scale Logical centring and scaling passed to
#'   [stats::prcomp()]. Both default to `TRUE`.
#' @param x_pc,y_pc Principal components to plot on the x and y axes.
#' @param group_colors Named character vector mapping groups to colours. `NULL`
#'   (default) uses the toolbox palette.
#' @param show_ellipse Whether to draw a per-group confidence ellipse.
#' @param ellipse_level Confidence level passed to [ggplot2::stat_ellipse()].
#' @param show_labels Whether to annotate each point with its sample label.
#' @param point_size Point size for samples.
#' @param point_alpha Point opacity from 0 to 1.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels. Defaults to `"PCn (xx.x%)"`.
#' @param font_family Font family used by the shared plot style. `"sans"` is
#'   the portable default.
#' @param style Optional object returned by [choose_plot_style()]. When
#'   supplied, it replaces the font and text-size defaults.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object with the PCA scores (`.PC1`, `.PC2`, `.Group`,
#'   `.Sample`) stored in `plot$data`.
#' @export
plot_pca <- function(
    data,
    sample_column = "Sample",
    group_column = "Group",
    expr_columns = NULL,
    center = TRUE,
    scale = TRUE,
    x_pc = 1,
    y_pc = 2,
    group_colors = NULL,
    show_ellipse = TRUE,
    ellipse_level = 0.95,
    show_labels = FALSE,
    point_size = 3,
    point_alpha = 0.8,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    font_family = "Times New Roman",
    style = NULL,
    output_file = NULL,
    figure_width = 8,
    figure_height = 7,
    dpi = 300) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(sample_column) || length(sample_column) != 1L ||
      is.na(sample_column) || !nzchar(sample_column)) {
    stop("`sample_column` must be one non-empty name.", call. = FALSE)
  }
  if (!is.character(group_column) || length(group_column) != 1L ||
      is.na(group_column) || !nzchar(group_column)) {
    stop("`group_column` must be one non-empty name.", call. = FALSE)
  }
  .assert_positive_number(point_size, "point_size")
  .assert_probability(point_alpha, "point_alpha")
  .assert_probability(ellipse_level, "ellipse_level")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  .assert_numeric_vector(x_pc, "x_pc", length_required = 1L, positive = TRUE)
  .assert_numeric_vector(y_pc, "y_pc", length_required = 1L, positive = TRUE)
  if (x_pc != floor(x_pc) || y_pc != floor(y_pc)) {
    stop("`x_pc` and `y_pc` must be whole numbers.", call. = FALSE)
  }
  x_pc <- as.integer(x_pc)
  y_pc <- as.integer(y_pc)
  .assert_flag(center, "center")
  .assert_flag(scale, "scale")
  .assert_flag(show_ellipse, "show_ellipse")
  .assert_flag(show_labels, "show_labels")

  missing_cols <- setdiff(c(sample_column, group_column), names(data))
  if (length(missing_cols)) {
    stop(
      "`data` is missing column(s): ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(expr_columns)) {
    expr_columns <- setdiff(names(data), c(sample_column, group_column))
  }
  if (!is.character(expr_columns) || !length(expr_columns) ||
      anyNA(expr_columns) || any(!nzchar(expr_columns))) {
    stop("`expr_columns` must be a non-empty character vector.", call. = FALSE)
  }
  missing_expr <- setdiff(expr_columns, names(data))
  if (length(missing_expr)) {
    stop(
      "`expr_columns` missing column(s): ", paste(missing_expr, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(expr_columns) < 2L) {
    stop("At least two expression columns are required for PCA.", call. = FALSE)
  }

  expr <- data[, expr_columns, drop = FALSE]
  for (column in expr_columns) {
    if (!is.numeric(expr[[column]])) {
      stop("Expression column `", column, "` is not numeric.", call. = FALSE)
    }
  }
  expr_matrix <- as.matrix(expr)
  suppressWarnings(storage.mode(expr_matrix) <- "double")

  valid_rows <- stats::complete.cases(expr_matrix) &
    apply(expr_matrix, 1L, function(row) all(is.finite(row)))
  if (any(!valid_rows)) {
    warning(
      sprintf("%d row(s) with non-finite expression values were excluded.", sum(!valid_rows)),
      call. = FALSE
    )
  }
  expr_matrix <- expr_matrix[valid_rows, , drop = FALSE]
  if (nrow(expr_matrix) < 2L) {
    stop("At least two complete rows are required for PCA.", call. = FALSE)
  }

  pca_result <- stats::prcomp(expr_matrix, center = center, scale. = scale)
  n_components <- length(pca_result$sdev)
  if (n_components < max(x_pc, y_pc)) {
    stop(
      sprintf("PCA produced only %d component(s); cannot plot PC%d vs PC%d.",
              n_components, x_pc, y_pc),
      call. = FALSE
    )
  }

  scores <- as.data.frame(pca_result$x)
  group <- data[[group_column]][valid_rows]
  if (is.factor(group)) group <- as.character(group)
  if (!is.character(group) || anyNA(group) || any(!nzchar(group))) {
    stop("`", group_column, "` cannot contain missing or empty values.", call. = FALSE)
  }
  scores[[".Group"]] <- factor(group, unique(group))
  scores[[".Sample"]] <- data[[sample_column]][valid_rows]
  scores[[".PC1"]] <- scores[[paste0("PC", x_pc)]]
  scores[[".PC2"]] <- scores[[paste0("PC", y_pc)]]

  group_levels <- levels(scores$.Group)
  if (is.null(group_colors)) {
    palette <- .toolbox_group_palette()
    if (length(group_levels) > length(palette)) {
      group_colors <- grDevices::hcl.colors(length(group_levels), "Dark 3")
      names(group_colors) <- group_levels
    } else {
      group_colors <- stats::setNames(palette[seq_along(group_levels)], group_levels)
    }
  }
  if (is.null(names(group_colors)) || anyNA(group_colors) ||
      any(!nzchar(group_colors))) {
    stop("`group_colors` must be a named character vector.", call. = FALSE)
  }
  missing_colors <- setdiff(group_levels, names(group_colors))
  if (length(missing_colors)) {
    stop(
      "`group_colors` is missing group(s): ",
      paste(missing_colors, collapse = ", "), call. = FALSE
    )
  }
  if (inherits(try(grDevices::col2rgb(group_colors), silent = TRUE), "try-error")) {
    stop("`group_colors` must contain valid R colors.", call. = FALSE)
  }

  var_exp <- 100 * pca_result$sdev^2 / sum(pca_result$sdev^2)
  if (is.null(x_label)) {
    x_label <- sprintf("PC%d (%.1f%%)", x_pc, var_exp[x_pc])
  }
  if (is.null(y_label)) {
    y_label <- sprintf("PC%d (%.1f%%)", y_pc, var_exp[y_pc])
  }

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
    scores,
    ggplot2::aes(x = .PC1, y = .PC2, color = .Group)
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(
      values = group_colors,
      breaks = group_levels,
      drop = FALSE,
      name = group_column
    ) +
    ggplot2::labs(x = x_label, y = y_label, title = title) +
    style$ggplot_theme

  if (show_ellipse) {
    plot <- plot + ggplot2::stat_ellipse(
      level = ellipse_level, linewidth = 0.6, alpha = 0.7
    )
  }
  if (show_labels) {
    plot <- plot + ggplot2::geom_text(
      ggplot2::aes(label = .Sample), vjust = -1, size = 3, show.legend = FALSE
    )
  }

  attr(plot, "x_pc") <- x_pc
  attr(plot, "y_pc") <- y_pc
  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

.pca_save_plot <- function(plot, output_file, width, height, dpi) {
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

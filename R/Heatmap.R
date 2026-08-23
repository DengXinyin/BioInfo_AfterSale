# Register common Microsoft font names for ComplexHeatmap's pre-device text
# measurements. Cairo/fontconfig still selects and embeds the real font glyphs.
.register_heatmap_font_metrics <- function(font_families) {
  metric_sources <- c("Times New Roman" = "Times", "Arial" = "Helvetica")
  registered <- names(grDevices::pdfFonts())
  for (font_family in unique(font_families)) {
    if (font_family %in% registered || !font_family %in% names(metric_sources)) next
    source_family <- unname(metric_sources[[font_family]])
    source_metrics <- grDevices::pdfFonts(source_family)[[1]]
    do.call(grDevices::pdfFonts, stats::setNames(list(source_metrics), font_family))
    registered <- c(registered, font_family)
  }
  invisible(NULL)
}

#' Draw a configurable ComplexHeatmap heatmap
#'
#' Draws an expression or Z-score matrix with optional sample-group
#' annotation. The Z-score legend and group legend are placed on the same side
#' of the heatmap; by default the Z-score legend is above the group legend on
#' the right. Row labels default to italic because they commonly represent gene
#' symbols.
#'
#' @param matrix Numeric matrix or data frame. Rows normally represent genes
#'   and columns represent samples.
#' @param group Optional sample-group vector. A named vector is reordered to
#'   match matrix column names; an unnamed vector must already be in column
#'   order.
#' @param group_colors Optional named character vector mapping groups to colors.
#' @param scale Scaling method: `"row"` for per-gene Z-scores, `"column"` for
#'   per-sample Z-scores, or `"none"` for an already transformed matrix.
#' @param show_row_names,show_column_names Whether to display row and column
#'   names.
#' @param cluster_rows,cluster_columns Whether to cluster rows and columns.
#' @param row_names_side,column_names_side Label sides accepted by
#'   [ComplexHeatmap::Heatmap()].
#' @param row_names_font_family,column_names_font_family Font families for row
#'   and column names.
#' @param row_names_font_size,column_names_font_size Font sizes in points.
#' @param row_names_italic Whether row names, usually gene symbols, are italic.
#' @param column_names_rot Numeric column-name rotation in degrees. Use `0` for
#'   horizontal labels and `90` for vertical labels; arbitrary angles are
#'   supported.
#' @param title Optional heatmap title.
#' @param title_position Title side: `"top"`, `"bottom"`, `"left"`, or
#'   `"right"`.
#' @param title_font_family,title_font_size,title_font_face Title font family,
#'   point size, and face (`"plain"`, `"bold"`, `"italic"`, or
#'   `"bold.italic"`).
#' @param zscore_legend_title,group_legend_title Legend titles.
#' @param show_zscore_legend,show_group_legend Whether to show each legend.
#' @param legend_side Side shared by both legends. The default is `"right"`.
#' @param group_annotation_side Place the column-group annotation above or
#'   below the heatmap.
#' @param zscore_breaks Three increasing numeric color breakpoints.
#' @param heatmap_colors Three colors corresponding to `zscore_breaks`.
#' @param na_color Color used for missing values.
#' @param border Whether to draw a border around the heatmap body.
#' @param output_file Optional PDF or PNG output filename.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#' @param draw_plot Whether to draw to the current graphics device when
#'   `output_file` is `NULL`. Set to `FALSE` to only return the Heatmap object.
#' @param ... Additional arguments passed to [ComplexHeatmap::Heatmap()].
#'
#' @return Invisibly returns the ComplexHeatmap object. When the plot is drawn,
#'   the returned object is the value from [ComplexHeatmap::draw()].
#' @export
Heatmap_plot <- function(
    matrix,
    group = NULL,
    group_colors = NULL,
    scale = c("row", "none", "column"),
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_side = c("right", "left"),
    column_names_side = c("bottom", "top"),
    row_names_font_family = "Times New Roman",
    column_names_font_family = "Times New Roman",
    row_names_font_size = 10,
    column_names_font_size = 12,
    row_names_italic = TRUE,
    column_names_rot = 0,
    title = NULL,
    title_position = c("top", "bottom", "left", "right"),
    title_font_family = "Times New Roman",
    title_font_size = 18,
    title_font_face = c("plain", "bold", "italic", "bold.italic"),
    zscore_legend_title = "Z-score",
    group_legend_title = "Group",
    show_zscore_legend = TRUE,
    show_group_legend = TRUE,
    legend_side = c("right", "left", "top", "bottom"),
    group_annotation_side = c("top", "bottom"),
    zscore_breaks = c(-2, 0, 2),
    heatmap_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
    na_color = "#BDBDBD",
    border = TRUE,
    output_file = NULL,
    figure_width = 8,
    figure_height = 8,
    dpi = 300,
    draw_plot = is.null(output_file),
    ...) {

  scale <- match.arg(scale)
  if (missing(zscore_legend_title) && scale == "none") {
    zscore_legend_title <- "Value"
  }
  row_names_side <- match.arg(row_names_side)
  column_names_side <- match.arg(column_names_side)
  title_position <- match.arg(title_position)
  title_font_face <- match.arg(title_font_face)
  legend_side <- match.arg(legend_side)
  group_annotation_side <- match.arg(group_annotation_side)

  for (argument in c(
    "show_row_names", "show_column_names", "cluster_rows", "cluster_columns",
    "row_names_italic", "show_zscore_legend", "show_group_legend", "border",
    "draw_plot"
  )) {
    .assert_flag(get(argument), argument)
  }
  .assert_positive_number(row_names_font_size, "row_names_font_size")
  .assert_positive_number(column_names_font_size, "column_names_font_size")
  .assert_positive_number(title_font_size, "title_font_size")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  .assert_numeric_vector(column_names_rot, "column_names_rot", length_required = 1L)
  .assert_numeric_vector(zscore_breaks, "zscore_breaks", length_required = 3L)
  if (is.unsorted(zscore_breaks, strictly = TRUE)) {
    stop("`zscore_breaks` must be strictly increasing.", call. = FALSE)
  }

  for (argument in c(
    "row_names_font_family", "column_names_font_family", "title_font_family",
    "zscore_legend_title", "group_legend_title", "na_color"
  )) {
    value <- get(argument)
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("`", argument, "` must be one non-empty character value.", call. = FALSE)
    }
  }
  if (!is.character(heatmap_colors) || length(heatmap_colors) != 3L ||
      anyNA(heatmap_colors) || any(!nzchar(heatmap_colors))) {
    stop("`heatmap_colors` must contain exactly three colors.", call. = FALSE)
  }
  if (inherits(try(grDevices::col2rgb(heatmap_colors), silent = TRUE), "try-error")) {
    stop("`heatmap_colors` must contain valid R colors.", call. = FALSE)
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("`title` must be NULL or one character value.", call. = FALSE)
  }
  .register_heatmap_font_metrics(c(
    row_names_font_family, column_names_font_family, title_font_family
  ))

  heat_matrix <- as.matrix(matrix)
  suppressWarnings(storage.mode(heat_matrix) <- "double")
  if (!length(heat_matrix) || nrow(heat_matrix) < 1L || ncol(heat_matrix) < 1L) {
    stop("`matrix` must contain at least one row and one column.", call. = FALSE)
  }
  if (all(!is.finite(heat_matrix))) {
    stop("`matrix` contains no finite numeric values.", call. = FALSE)
  }

  if (scale == "row") {
    row_sd <- apply(heat_matrix, 1L, stats::sd, na.rm = TRUE)
    invalid_rows <- !is.finite(row_sd) | row_sd == 0
    if (any(invalid_rows)) {
      removed_rows <- rownames(heat_matrix)[invalid_rows]
      removed_rows <- removed_rows[!is.na(removed_rows)]
      warning(
        sprintf(
          "%d zero-variance row(s) were removed before row Z-score scaling%s.",
          sum(invalid_rows),
          if (length(removed_rows)) paste0(": ", paste(removed_rows, collapse = ", ")) else ""
        ),
        call. = FALSE
      )
      heat_matrix <- heat_matrix[!invalid_rows, , drop = FALSE]
      if (!nrow(heat_matrix)) {
        stop("No variable rows remain for row Z-score scaling.", call. = FALSE)
      }
    }
    heat_matrix <- t(base::scale(t(heat_matrix)))
  } else if (scale == "column") {
    heat_matrix <- base::scale(heat_matrix)
  }
  heat_matrix[is.nan(heat_matrix)] <- 0

  color_function <- circlize::colorRamp2(zscore_breaks, heatmap_colors)
  manual_legends <- list()
  if (show_zscore_legend) {
    manual_legends[[length(manual_legends) + 1L]] <- ComplexHeatmap::Legend(
      title = zscore_legend_title,
      col_fun = color_function,
      at = zscore_breaks,
      title_position = "topcenter"
    )
  }

  annotation <- NULL
  if (!is.null(group)) {
    if (!is.atomic(group)) stop("`group` must be an atomic vector.", call. = FALSE)
    if (!is.null(names(group))) {
      if (is.null(colnames(heat_matrix))) {
        stop("Named `group` requires matrix column names.", call. = FALSE)
      }
      missing_groups <- setdiff(colnames(heat_matrix), names(group))
      if (length(missing_groups)) {
        stop(
          "Named `group` is missing sample(s): ",
          paste(missing_groups, collapse = ", "), call. = FALSE
        )
      }
      group <- group[colnames(heat_matrix)]
    } else if (length(group) != ncol(heat_matrix)) {
      stop("Unnamed `group` must have one value per matrix column.", call. = FALSE)
    }
    if (anyNA(group) || any(!nzchar(as.character(group)))) {
      stop("`group` cannot contain missing or empty values.", call. = FALSE)
    }
    group <- if (is.factor(group)) group else factor(as.character(group), unique(group))
    group_levels <- levels(droplevels(group))
    if (is.null(group_colors)) {
      group_colors <- grDevices::hcl.colors(length(group_levels), "Dark 3")
      names(group_colors) <- group_levels
    }
    if (!is.character(group_colors) || is.null(names(group_colors)) ||
        anyNA(group_colors) || any(!nzchar(group_colors))) {
      stop("`group_colors` must be a named character vector.", call. = FALSE)
    }
    missing_colors <- setdiff(group_levels, names(group_colors))
    if (length(missing_colors)) {
      stop(
        "`group_colors` is missing group(s): ",
        paste(missing_colors, collapse = ", "), call. = FALSE
      )
    }
    annotation <- ComplexHeatmap::HeatmapAnnotation(
      Group = group,
      col = list(Group = group_colors[group_levels]),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    if (show_group_legend) {
      manual_legends[[length(manual_legends) + 1L]] <- ComplexHeatmap::Legend(
        title = group_legend_title,
        labels = group_levels,
        legend_gp = grid::gpar(fill = unname(group_colors[group_levels])),
        direction = "vertical",
        title_position = "topleft"
      )
    }
  }

  title_args <- list(
    column_title = NULL,
    row_title = NULL,
    column_title_gp = grid::gpar(
      fontfamily = title_font_family, fontsize = title_font_size,
      fontface = title_font_face
    ),
    row_title_gp = grid::gpar(
      fontfamily = title_font_family, fontsize = title_font_size,
      fontface = title_font_face
    )
  )
  if (!is.null(title) && nzchar(title)) {
    if (title_position %in% c("top", "bottom")) {
      title_args$column_title <- title
      title_args$column_title_side <- title_position
    } else {
      title_args$row_title <- title
      title_args$row_title_side <- title_position
    }
  }

  annotation_args <- if (group_annotation_side == "top") {
    list(top_annotation = annotation)
  } else {
    list(bottom_annotation = annotation)
  }
  heatmap_args <- c(
    list(
      matrix = heat_matrix,
      name = zscore_legend_title,
      col = color_function,
      na_col = na_color,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_side = row_names_side,
      column_names_side = column_names_side,
      row_names_gp = grid::gpar(
        fontfamily = row_names_font_family,
        fontsize = row_names_font_size,
        fontface = if (row_names_italic) "italic" else "plain"
      ),
      column_names_gp = grid::gpar(
        fontfamily = column_names_font_family,
        fontsize = column_names_font_size
      ),
      column_names_rot = column_names_rot,
      show_heatmap_legend = FALSE,
      heatmap_legend_param = list(
        title = zscore_legend_title,
        title_position = "topcenter",
        at = zscore_breaks
      ),
      border = border
    ),
    annotation_args,
    title_args,
    list(...)
  )
  heatmap <- do.call(ComplexHeatmap::Heatmap, heatmap_args)
  draw_args <- list(
    object = heatmap,
    heatmap_legend_side = legend_side,
    annotation_legend_side = legend_side,
    merge_legends = TRUE,
    legend_grouping = "original",
    heatmap_legend_list = manual_legends,
    annotation_legend_list = list()
  )
  attr(heatmap, "aftersale_draw_args") <- draw_args[-1]

  draw_heatmap <- function() do.call(ComplexHeatmap::draw, draw_args)
  drawn <- NULL
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1L ||
        is.na(output_file) || !nzchar(output_file)) {
      stop("`output_file` must be NULL or one non-empty filename.", call. = FALSE)
    }
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      stop("Output directory does not exist: ", output_dir, call. = FALSE)
    }
    extension <- tolower(tools::file_ext(output_file))
    if (!extension %in% c("pdf", "png")) {
      stop("`output_file` must end in .pdf or .png.", call. = FALSE)
    }
    if (extension == "pdf") {
      grDevices::cairo_pdf(
        output_file, width = figure_width, height = figure_height,
        family = title_font_family, onefile = FALSE
      )
    } else {
      grDevices::png(
        output_file, width = figure_width, height = figure_height,
        units = "in", res = dpi, type = "cairo", bg = "white"
      )
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    drawn <- draw_heatmap()
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
  } else if (draw_plot) {
    drawn <- draw_heatmap()
  }

  invisible(if (is.null(drawn)) heatmap else drawn)
}

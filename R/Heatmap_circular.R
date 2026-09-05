# Circular (ring) clustered heatmap for an expression matrix.
#
# Reference: "如何使用DeepSeek绘制好看的环形热图" (circos.heatmap pipeline,
#   gene x sample expression matrix, row Z-score scaling, inner dendrogram,
#   outer gene labels, legend, adjustable start/gap angles).
#
# Layout notes verified against circlize internals:
#   * circos.heatmap() places matrix ROWS (genes) along the circle and matrix
#     COLUMNS (samples) as radial layers inside one track.
#   * The dendrogram is drawn by circos.heatmap() itself; choose
#     dend.side = "inside" (default) or "outside".
#   * Sample (column) names are not drawn by circos.heatmap(); we add a
#     dedicated radial text layer on the heatmap track with circos.text().

#' Plot a circular (ring) clustered heatmap
#'
#' Draws a gene-by-sample expression matrix as a circular heatmap with
#' per-row Z-score scaling (optional), hierarchical clustering of genes, an
#' inner (or outer) dendrogram, gene labels on the outer rim, sample labels
#' arranged radially, and a continuous colour legend. The implementation is
#' based on `circlize::circos.heatmap()`; the adjustable start angle and gap
#' angle follow the shared circular-heatmap tutorial workflow.
#'
#' @param title Optional heatmap title.
#' @param style A style object created by [choose_plot_style()].
#' @param matrix Numeric matrix or data frame of expression values. Rows are
#'   genes (drawn along the circle), columns are samples (drawn as radial
#'   layers).
#' @param scale Row/column scaling before plotting: `"row"` applies per-gene
#'   Z-scores (default, recommended for expression matrices), `"column"`
#'   applies per-sample Z-scores, and `"none"` plots values as supplied.
#' @param cluster Whether to cluster genes (rows). When `TRUE`, an average
#'   linkage dendrogram on Euclidean distances is shown inside (or outside,
#'   see `dend_side`) the ring.
#' @param dend_side Side for the gene dendrogram: `"inside"` or `"outside"`.
#' @param rownames_side Side for gene labels: `"inside"` or `"outside"`.
#'   Defaults to the opposite of `dend_side`.
#' @param show_rownames Whether to draw gene (row) labels.
#' @param rownames_cex Gene-label size relative to the default device size.
#' @param rownames_col Gene-label colour.
#' @param show_colnames Whether to draw sample (column) labels as a radial
#'   text layer inside the heatmap track.
#' @param colnames_cex Sample-label size relative to the default device size.
#' @param start_degree Degree position of the first sector (ring start).
#' @param gap_degree Angular gap between adjacent sectors.
#' @param track_height Relative height of the heatmap track.
#' @param dend_track_height Relative height of the dendrogram track.
#' @param cell_border Border colour of every heatmap cell.
#' @param cell_lwd Border line width of every heatmap cell.
#' @param zscore_breaks Three increasing numeric colour breakpoints used with
#'   `heatmap_colors`.
#' @param heatmap_colors Three colours corresponding to `zscore_breaks`
#'   (low / mid / high). When `scale = "none"` and these are not supplied, the
#'   breaks are taken from the observed range of the matrix.
#' @param na_color Colour for missing values.
#' @param legend_title Title of the continuous colour legend.
#' @param show_legend Whether to draw the continuous colour legend.
#' @param output_file Optional PDF or PNG output filename. For a circular
#'   heatmap, PDF is the recommended format (vector output).
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#' @param draw_plot Whether to draw to the current graphics device when
#'   `output_file` is `NULL`.
#' @param ... Additional arguments passed to [circlize::circos.heatmap()].
#'
#' @return Invisibly returns the scaled matrix that was plotted (rows in
#'   dendrogram order when `cluster = TRUE`).
#' @export
plot_circular_heatmap <- function(
    matrix,
    scale = c("row", "none", "column"),
    cluster = TRUE,
    dend_side = c("inside", "outside"),
    rownames_side = NULL,
    show_rownames = TRUE,
    rownames_cex = 0.7,
    rownames_col = "black",
    show_colnames = TRUE,
    colnames_cex = 0.7,
    start_degree = 45,
    gap_degree = 45,
    track_height = 0.3,
    dend_track_height = 0.15,
    cell_border = "white",
    cell_lwd = 0.7,
    zscore_breaks = c(-2, 0, 2),
    heatmap_colors = c("#A5CC26", "white", "#FF7BAC"),
    na_color = "#BDBDBD",
    legend_title = "Exp",
    show_legend = TRUE,
    title = NULL,
    style = NULL,
    output_file = NULL,
    figure_width = 10,
    figure_height = 10,
    dpi = 300,
    draw_plot = is.null(output_file),
    ...) {

  scale <- match.arg(scale)
  dend_side <- match.arg(dend_side)
  rownames_side <- rownames_side %||%
    (if (dend_side == "inside") "outside" else "inside")
  rownames_side <- match.arg(rownames_side, c("inside", "outside"))
  if (identical(dend_side, rownames_side)) {
    stop("`dend_side` and `rownames_side` must differ.", call. = FALSE)
  }

  for (argument in c(
    "show_rownames", "show_colnames", "cluster", "show_legend", "draw_plot"
  )) {
    .assert_flag(get(argument), argument)
  }
  .assert_positive_number(rownames_cex, "rownames_cex")
  .assert_positive_number(colnames_cex, "colnames_cex")
  .assert_positive_number(start_degree, "start_degree")
  .assert_positive_number(gap_degree, "gap_degree")
  .assert_positive_number(track_height, "track_height")
  .assert_positive_number(dend_track_height, "dend_track_height")
  .assert_positive_number(cell_lwd, "cell_lwd")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  for (argument in c(
    "rownames_col", "cell_border", "na_color", "legend_title"
  )) {
    value <- get(argument)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("`", argument, "` must be one non-empty character value.",
           call. = FALSE)
    }
  }
  if (!is.character(heatmap_colors) || length(heatmap_colors) != 3L ||
      anyNA(heatmap_colors) || any(!nzchar(heatmap_colors))) {
    stop("`heatmap_colors` must contain exactly three colors.", call. = FALSE)
  }
  if (inherits(try(grDevices::col2rgb(heatmap_colors), silent = TRUE),
               "try-error")) {
    stop("`heatmap_colors` must contain valid R colors.", call. = FALSE)
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("`title` must be NULL or one character value.", call. = FALSE)
  }

  style <- .plot_style_or_default(style)
  .register_heatmap_font_metrics(style$global$font_family)

  heat_matrix <- as.matrix(matrix)
  suppressWarnings(storage.mode(heat_matrix) <- "double")
  if (!length(heat_matrix) || nrow(heat_matrix) < 1L ||
      ncol(heat_matrix) < 1L) {
    stop("`matrix` must contain at least one row and one column.",
         call. = FALSE)
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
          if (length(removed_rows)) {
            paste0(": ", paste(removed_rows, collapse = ", "))
          } else ""
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

  observed <- range(heat_matrix[is.finite(heat_matrix)], na.rm = TRUE)
  if (diff(observed) == 0) {
    stop("`matrix` has a constant value range; cannot build a colour scale.",
         call. = FALSE)
  }
  breaks <- zscore_breaks
  if (scale == "none" && missing(heatmap_colors)) {
    breaks <- seq(observed[1], observed[2], length.out = 3L)
    heatmap_colors <- c("#2166AC", "#F7F7F7", "#B2182B")
  }
  if (is.unsorted(breaks, strictly = TRUE)) {
    stop("`zscore_breaks` must be strictly increasing.", call. = FALSE)
  }
  col_fun <- circlize::colorRamp2(breaks, heatmap_colors)

  # n-sized font metric helper: recycles a single value like circlize does.
  n_genes <- nrow(heat_matrix)
  n_samples <- ncol(heat_matrix)

  draw_one <- function() {
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(
      start.degree = start_degree,
      gap.degree = gap_degree,
      points.overflow.warning = FALSE
    )
    circlize::circos.heatmap(
      heat_matrix,
      col = col_fun,
      na.col = na_color,
      cluster = cluster,
      clustering.method = "average",
      distance.method = "euclidean",
      dend.side = dend_side,
      dend.track.height = dend_track_height,
      rownames.side = if (show_rownames) rownames_side else "none",
      rownames.cex = rownames_cex,
      rownames.font = 1,
      rownames.col = rownames_col,
      bg.border = NA,
      cell.border = if (is.na(cell_border)) NA else cell_border,
      cell.lwd = cell_lwd,
      track.height = track_height,
      ...
    )

    if (show_colnames) {
      heat_track <- 2L
      if (!show_rownames || rownames_side != "outside") heat_track <- 1L
      circlize::circos.track(
        track.index = heat_track,
        panel.fun = function(x, y) {
          if (circlize::get.cell.meta.data("sector.numeric.index") == 1L) {
            cn <- colnames(heat_matrix)
            n <- length(cn)
            circlize::circos.text(
              rep(circlize::get.cell.meta.data("cell.xlim")[2], n) +
                circlize::convert_x(2, "mm"),
              1:n - 0.5, cn,
              cex = colnames_cex, adj = c(0, 0.6), facing = "inside",
              font = 1
            )
          }
        },
        bg.border = NA
      )
    }

    if (!is.null(title) && nzchar(title)) {
      graphics::title(
        main = title,
        line = -3,
        cex.main = style$text$title$size / 12,
        font.main = if (style$text$title$bold) 2 else 1,
        family = style$global$font_family
      )
    }
    if (show_legend) {
      lg <- ComplexHeatmap::Legend(
        title = legend_title,
        col_fun = col_fun,
        at = breaks,
        labels = format(breaks, trim = TRUE),
        direction = "vertical",
        title_position = "topcenter",
        title_gp = grid::gpar(
          fontsize = style$text$legend_title$size,
          fontfamily = style$global$font_family
        ),
        labels_gp = grid::gpar(
          fontsize = style$text$legend_text$size,
          fontfamily = style$global$font_family
        )
      )
      ComplexHeatmap::draw(
        lg,
        x = grid::unit(0.92, "npc"),
        y = grid::unit(0.5, "npc"),
        just = c("right", "top")
      )
    }
  }

  drawn <- NULL
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1L ||
        is.na(output_file) || !nzchar(output_file)) {
      stop("`output_file` must be NULL or one non-empty filename.",
           call. = FALSE)
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
        family = style$global$font_family, onefile = FALSE
      )
    } else {
      grDevices::png(
        output_file, width = figure_width, height = figure_height,
        units = "in", res = dpi, type = "cairo", bg = "white"
      )
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_one()
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
  } else if (draw_plot) {
    draw_one()
  }

  invisible(heat_matrix)
}

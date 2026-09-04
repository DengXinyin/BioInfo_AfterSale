#' Plot a circlize circular cross-correlation heatmap
#'
#' Rows of `x` and `y` are matched by sample row names. The correlation
#' matrix is rendered as concentric rings with taxa/features as sectors,
#' average-linkage Euclidean clustering, an inner dendrogram, and a shared
#' continuous correlation legend.
#'
#' @inheritParams plot_dot_heatmap
#' @param cluster Whether to cluster row sectors with average-linkage
#'   Euclidean clustering.
#' @return Invisibly returns the correlation matrix.
#' @export
plot_circular_heatmap <- function(x, y, method = "spearman", title = NULL,
                                  style = NULL, output_file = NULL, width = NULL,
                                  height = NULL, dpi = NULL, cluster = TRUE) {
  style <- .plot_style_or_default(style)
  if (!is.logical(cluster) || length(cluster) != 1L || is.na(cluster)) stop("`cluster` must be TRUE or FALSE.", call. = FALSE)
  rho <- .correlation_matrix(x, y, method)
  rho <- rho[rowSums(is.finite(rho)) > 0, , drop = FALSE]
  if (!nrow(rho) || !ncol(rho)) stop("The correlation matrix contains no finite values.", call. = FALSE)
  if (!is.null(output_file) && (!is.character(output_file) || length(output_file) != 1L || is.na(output_file) || !nzchar(output_file))) stop("`output_file` must be NULL or one filename.", call. = FALSE)
  width <- width %||% style$global$figure_width; height <- height %||% style$global$figure_height; dpi <- dpi %||% style$global$dpi
  if (!is.null(output_file) && !dir.exists(dirname(output_file))) stop("Output directory does not exist: ", dirname(output_file), call. = FALSE)
  draw_one <- function() {
    if (!requireNamespace("circlize", quietly = TRUE) || !requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Packages `circlize` and `ComplexHeatmap` are required for the circular heatmap.", call. = FALSE)
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
    equal_level_dend <- function(dend, ...) {
      set_levels <- function(node) {
        if (is.leaf(node)) { attr(node, "height") <- 0; return(node) }
        for (i in seq_along(node)) node[[i]] <- set_levels(node[[i]])
        attr(node, "height") <- max(vapply(node, function(a) attr(a, "height"), numeric(1))) + 1
        node
      }
      set_levels(dend)
    }
    n <- nrow(rho); cex <- max(.55, min(1.15, 19.4 / max(n, 22) + 6 / 12))
    circlize::circos.par(start.degree = 90, gap.after = 75, cell.padding = c(0, 0, 0, 0), track.margin = c(.001, .001), canvas.xlim = c(-2.5, 2.5), canvas.ylim = c(-2.5, 2.5), points.overflow.warning = FALSE)
    circlize::circos.heatmap(rho, col = col_fun, na.col = "white", cluster = cluster, clustering.method = "average", distance.method = "euclidean", dend.callback = equal_level_dend, dend.side = "inside", dend.track.height = .12, rownames.side = "outside", rownames.cex = cex, rownames.font = 1, rownames.col = "black", bg.border = NA, cell.border = "white", cell.lwd = .25, track.height = .60)
    if (!is.null(title)) graphics::title(main = title, line = -3, cex.main = style$text$title$size / 12, font.main = if (style$text$title$bold) 2 else 1, family = style$global$font_family)
    lg <- ComplexHeatmap::Legend(title = paste0(method, " r"), col_fun = col_fun, at = c(-1, -.5, 0, .5, 1), labels = c("-1", "-0.5", "0", "0.5", "1"), direction = "vertical", title_position = "topcenter", title_gp = grid::gpar(fontsize = style$text$legend_title$size), labels_gp = grid::gpar(fontsize = style$text$legend_text$size))
    ComplexHeatmap::draw(lg, x = grid::unit(.98, "npc"), y = grid::unit(.98, "npc"), just = c("right", "top"))
  }
  if (!is.null(output_file)) {
    ext <- tolower(tools::file_ext(output_file)); if (ext != "pdf") stop("`output_file` must end in .pdf for a circular heatmap.", call. = FALSE)
    grDevices::cairo_pdf(output_file, width = width, height = height, family = style$global$font_family); on.exit(grDevices::dev.off(), add = TRUE); draw_one(); grDevices::dev.off(); on.exit(NULL, add = FALSE)
  } else draw_one()
  invisible(rho)
}

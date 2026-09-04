#' Draw a protein interaction or co-expression network
#'
#' @param edges Edge table.
#' @param nodes Node table.
#' @param source,target Edge endpoint column names.
#' @param node_id,status Node ID and category column names.
#' @param weight Optional numeric edge-weight column.
#' @param colors Optional named node-category colors.
#' @param seed Random seed for the force-directed layout.
#' @param label_top Number of highest-degree node labels displayed.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
ppi_network <- function(edges, nodes, source = "Node1", target = "Node2",
                        node_id = "Gene", status = "Status", weight = NULL,
                        colors = c(Up = "#E15759", Down = "#59A14F", NoSig = "#A0A0A0"),
                        seed = 2026, label_top = 0, title = "Interaction network",
                        style = NULL, output_file = NULL, width = NULL,
                        height = NULL, dpi = NULL) {
  .sanshu_require_columns(edges, c(source, target, weight))
  .sanshu_require_columns(nodes, c(node_id, status))
  if (!is.null(weight)) .sanshu_numeric(edges, weight)
  if (anyDuplicated(nodes[[node_id]])) stop("Node IDs must be unique.", call. = FALSE)
  endpoints <- unique(c(as.character(edges[[source]]), as.character(edges[[target]])))
  missing <- setdiff(endpoints, as.character(nodes[[node_id]]))
  if (length(missing)) {
    stop("`nodes` is missing edge endpoint(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  .assert_positive_number(seed, "seed")
  if (!is.numeric(label_top) || length(label_top) != 1L || is.na(label_top) || label_top < 0) {
    stop("`label_top` must be a non-negative number.", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network layout.", call. = FALSE)
  }
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  edge_data <- data.frame(
    From = as.character(edges[[source]]), To = as.character(edges[[target]]),
    Weight = if (is.null(weight)) 1 else edges[[weight]], stringsAsFactors = FALSE
  )
  node_data <- data.frame(
    Gene = as.character(nodes[[node_id]]), Status = as.character(nodes[[status]]),
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(edge_data[, c("From", "To")],
                                          directed = FALSE, vertices = node_data$Gene)
  set.seed(as.integer(seed))
  coordinates <- igraph::layout_with_fr(graph, weights = edge_data$Weight)
  node_data$X <- coordinates[, 1]
  node_data$Y <- coordinates[, 2]
  node_data$Degree <- igraph::degree(graph, v = node_data$Gene)
  names(node_data$X) <- names(node_data$Y) <- node_data$Gene
  edge_data$X <- node_data$X[match(edge_data$From, node_data$Gene)]
  edge_data$Y <- node_data$Y[match(edge_data$From, node_data$Gene)]
  edge_data$Xend <- node_data$X[match(edge_data$To, node_data$Gene)]
  edge_data$Yend <- node_data$Y[match(edge_data$To, node_data$Gene)]
  statuses <- unique(node_data$Status)
  palette <- .sanshu_palette(statuses, colors, z$style)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_data, ggplot2::aes(
      x = .data$X, y = .data$Y, xend = .data$Xend, yend = .data$Yend,
      linewidth = .data$Weight
    ), color = "#BDBDBD", alpha = 0.2, show.legend = FALSE) +
    ggplot2::scale_linewidth(range = c(0.1, 0.8)) +
    ggplot2::geom_point(data = node_data, ggplot2::aes(
      x = .data$X, y = .data$Y, color = .data$Status,
      size = sqrt(.data$Degree + 1)
    ), alpha = 0.9) +
    ggplot2::scale_color_manual(values = palette, breaks = statuses, name = NULL) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
  if (label_top > 0) {
    labelled <- utils::head(node_data[order(node_data$Degree, decreasing = TRUE), , drop = FALSE],
                            as.integer(label_top))
    p <- p + ggplot2::geom_text(
      data = labelled, ggplot2::aes(.data$X, .data$Y, label = .data$Gene),
      family = z$style$global$font_family,
      size = z$style$text$data_label$size / 3.2, vjust = -0.8
    )
  }
  .sanshu_finish(p, z)
}

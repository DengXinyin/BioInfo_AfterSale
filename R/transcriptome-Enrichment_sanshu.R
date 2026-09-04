#' Draw an enrichment butterfly chart
#'
#' @param data Enrichment table.
#' @param term,up,down Column names for term and regulated-gene counts.
#' @param colors Named colors for `Up` and `Down`.
#' @param top_n Number of highest-total terms shown.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
enrichment_butterfly <- function(data, term = "Term", up = "Up", down = "Down",
                                 colors = c(Up = "#D5695D", Down = "#65A479"),
                                 top_n = 30, title = NULL, style = NULL,
                                 output_file = NULL, width = NULL, height = NULL,
                                 dpi = NULL) {
  .sanshu_require_columns(data, c(term, up, down))
  .sanshu_numeric(data, c(up, down))
  if (any(data[[up]] < 0) || any(data[[down]] < 0)) {
    stop("Up/down counts cannot be negative.", call. = FALSE)
  }
  .assert_positive_number(top_n, "top_n")
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- data
  d$.total <- d[[up]] + d[[down]]
  d <- utils::head(d[order(d$.total, decreasing = TRUE), , drop = FALSE], as.integer(top_n))
  terms <- rev(as.character(d[[term]]))
  long <- rbind(
    data.frame(Term = as.character(d[[term]]), Direction = "Up", Count = d[[up]]),
    data.frame(Term = as.character(d[[term]]), Direction = "Down", Count = -d[[down]])
  )
  long$Term <- factor(long$Term, levels = terms)
  palette <- .sanshu_palette(c("Up", "Down"), colors, z$style)
  maximum <- max(abs(long$Count))
  breaks <- pretty(c(-maximum, maximum), n = 5)
  p <- ggplot2::ggplot(long, ggplot2::aes(
    x = .data$Term, y = .data$Count, fill = .data$Direction
  )) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(breaks = breaks, labels = abs(breaks)) +
    ggplot2::scale_fill_manual(values = palette, breaks = c("Up", "Down"), name = NULL) +
    ggplot2::labs(x = NULL, y = "Mapping numbers", title = title) +
    z$style$ggplot_theme
  .sanshu_finish(p, z)
}

#' Draw a multi-track enrichment circle
#'
#' The tracks reproduce the source layout from outside to inside: term ID and
#' category band, gene count colored by `-log10(Padj)`, proportional up/down
#' sectors, and rich-factor bars.
#'
#' @param data Enrichment table.
#' @param id,category,count,up,down,rich_factor,padj Column names.
#' @param category_colors Named category colors.
#' @param top_n Number of lowest-Padj terms shown.
#' @inheritParams enrichment_butterfly
#' @return A ggplot object.
#' @export
enrichment_circle <- function(
    data, id = "ID", category = "Category", count = "Count", up = "Up",
    down = "Down", rich_factor = "RichFactor", padj = "Padj",
    category_colors = NULL, top_n = 15, title = NULL, style = NULL,
    output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  .sanshu_require_columns(data, c(id, category, count, up, down, rich_factor, padj))
  .sanshu_numeric(data, c(count, up, down, rich_factor, padj))
  if (any(data[[count]] < 0 | data[[up]] < 0 | data[[down]] < 0 |
          data[[rich_factor]] < 0 | data[[padj]] <= 0 | data[[padj]] > 1)) {
    stop("Counts/rich factors must be non-negative and Padj must satisfy 0 < Padj <= 1.", call. = FALSE)
  }
  .assert_positive_number(top_n, "top_n")
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  d <- utils::head(data[order(data[[padj]], decreasing = FALSE), , drop = FALSE], as.integer(top_n))
  if (anyDuplicated(d[[id]])) stop("Displayed enrichment IDs must be unique.", call. = FALSE)
  d$.x <- seq_len(nrow(d))
  d$.logp <- -log10(d[[padj]])
  d$.count_height <- 1.6 * d[[count]] / max(d[[count]], 1)
  d$.rich_height <- 1.5 * d[[rich_factor]] / max(d[[rich_factor]], 1e-12)
  categories <- unique(as.character(d[[category]]))
  palette <- .sanshu_palette(categories, category_colors, z$style)
  d$.category_color <- unname(palette[as.character(d[[category]])])
  sector_width <- 0.92
  total <- d[[up]] + d[[down]]
  down_fraction <- ifelse(total > 0, d[[down]] / total, 0.5)
  d$.down_xmin <- d$.x - sector_width / 2
  d$.down_xmax <- d$.down_xmin + sector_width * down_fraction
  d$.up_xmin <- d$.down_xmax
  d$.up_xmax <- d$.x + sector_width / 2
  d$.down_x <- (d$.down_xmin + d$.down_xmax) / 2
  d$.up_x <- (d$.up_xmin + d$.up_xmax) / 2
  angle <- 90 - 360 * (d$.x - 0.5) / nrow(d)
  d$.angle <- ifelse(angle < -90, angle + 180, angle)
  d$.hjust <- ifelse(angle < -90, 1, 0)
  legend_levels <- c(categories, "Up-regulated", "Down-regulated")
  legend_colors <- c(palette, "Up-regulated" = "#CC0000", "Down-regulated" = "#2457C5")
  legend_data <- data.frame(X = 1, Y = 0, Legend = factor(legend_levels, levels = legend_levels))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = d, ggplot2::aes(
      xmin = .data$.x - sector_width / 2, xmax = .data$.x + sector_width / 2,
      ymin = 9.9, ymax = 10.45
    ), fill = d$.category_color, color = "white", linewidth = 0.25) +
    ggplot2::geom_rect(data = d, ggplot2::aes(
      xmin = .data$.x - sector_width / 2, xmax = .data$.x + sector_width / 2,
      ymin = 8, ymax = 8 + .data$.count_height, fill = .data$.logp
    ), color = "white", linewidth = 0.25) +
    ggplot2::geom_text(data = d, ggplot2::aes(
      x = .data$.x, y = 8 + .data$.count_height / 2, label = .data[[count]],
      angle = .data$.angle
    ), color = "white", family = z$style$global$font_family,
    size = z$style$text$data_label$size / 4) +
    ggplot2::geom_rect(data = d, ggplot2::aes(
      xmin = .data$.down_xmin, xmax = .data$.down_xmax, ymin = 5.8, ymax = 6.55
    ), fill = "#2457C5", color = "white", linewidth = 0.2) +
    ggplot2::geom_rect(data = d, ggplot2::aes(
      xmin = .data$.up_xmin, xmax = .data$.up_xmax, ymin = 5.8, ymax = 6.55
    ), fill = "#CC0000", color = "white", linewidth = 0.2) +
    ggplot2::geom_text(data = d[d[[down]] > 0, , drop = FALSE], ggplot2::aes(
      x = .data$.down_x, y = 5.63, label = .data[[down]], angle = .data$.angle
    ), family = z$style$global$font_family,
    size = z$style$text$data_label$size / 4) +
    ggplot2::geom_text(data = d[d[[up]] > 0, , drop = FALSE], ggplot2::aes(
      x = .data$.up_x, y = 5.63, label = .data[[up]], angle = .data$.angle
    ), family = z$style$global$font_family,
    size = z$style$text$data_label$size / 4) +
    ggplot2::geom_rect(data = d, ggplot2::aes(
      xmin = .data$.x - sector_width / 2, xmax = .data$.x + sector_width / 2,
      ymin = 4, ymax = 4 + .data$.rich_height
    ), fill = d$.category_color, color = "white", linewidth = 0.2) +
    ggplot2::geom_text(data = d, ggplot2::aes(
      x = .data$.x, y = 10.62, label = .data[[id]], angle = .data$.angle,
      hjust = .data$.hjust
    ), family = z$style$global$font_family,
    size = z$style$text$axis_text$size / 3.2) +
    ggplot2::geom_point(data = legend_data, ggplot2::aes(
      x = .data$X, y = .data$Y, color = .data$Legend
    ), alpha = 0) +
    ggplot2::scale_fill_gradient(
      low = "#FF906F", high = "#861D30", name = "-log10(Padj)",
      breaks = .default_numeric_breaks(d$.logp, 3)
    ) +
    ggplot2::scale_color_manual(values = legend_colors, name = NULL, drop = FALSE) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1, override.aes = list(alpha = 1, shape = 15, size = 4)
      ),
      fill = ggplot2::guide_colorbar(
        order = 2, direction = "horizontal", title.position = "top",
        label.position = "bottom", barwidth = 4, barheight = 0.45
      )
    ) +
    ggplot2::coord_polar(theta = "x", clip = "off") +
    ggplot2::scale_x_continuous(limits = c(0.5, nrow(d) + 0.5)) +
    ggplot2::scale_y_continuous(limits = c(0, 11.5)) +
    ggplot2::labs(title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = c(0.5, 0.5),
                   legend.justification = c(0.5, 0.5),
                   legend.background = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_blank())
  attr(p, "category_colors") <- palette
  .sanshu_finish(p, z)
}

#' Draw a GO directed-acyclic graph from explicit nodes and edges
#'
#' @param nodes Node table.
#' @param edges Edge table directed from parent to child.
#' @param id,term,pvalue,gene_ratio Node column names.
#' @param parent,child Edge column names.
#' @param title Plot title.
#' @inheritParams enrichment_butterfly
#' @return A ggplot object.
#' @export
go_dag <- function(nodes, edges, id = "ID", term = "Term", pvalue = "Pvalue",
                   gene_ratio = "GeneRatio", parent = "Parent", child = "Child",
                   title = "GO enrichment DAG", style = NULL,
                   output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  .sanshu_require_columns(nodes, c(id, term, pvalue, gene_ratio))
  .sanshu_require_columns(edges, c(parent, child))
  .sanshu_numeric(nodes, pvalue)
  if (any(nodes[[pvalue]] <= 0 | nodes[[pvalue]] > 1)) {
    stop("Node P values must satisfy 0 < Pvalue <= 1.", call. = FALSE)
  }
  if (anyDuplicated(nodes[[id]])) stop("Node IDs must be unique.", call. = FALSE)
  node_ids <- as.character(nodes[[id]])
  endpoints <- unique(c(as.character(edges[[parent]]), as.character(edges[[child]])))
  if (!setequal(node_ids, endpoints)) {
    stop("Node IDs and edge endpoints must match exactly.", call. = FALSE)
  }
  z <- .sanshu_style_output(style, output_file, width, height, dpi)
  parents <- as.character(edges[[parent]])
  children <- as.character(edges[[child]])
  indegree <- table(factor(children, levels = node_ids))
  queue <- node_ids[indegree == 0]
  if (!length(queue)) stop("The graph has no root and may contain a cycle.", call. = FALSE)
  level <- stats::setNames(rep(0L, length(node_ids)), node_ids)
  visited <- character()
  while (length(queue)) {
    current <- queue[1]
    queue <- queue[-1]
    visited <- c(visited, current)
    targets <- children[parents == current]
    for (target in targets) {
      level[target] <- max(level[target], level[current] + 1L)
      indegree[target] <- indegree[target] - 1L
      if (indegree[target] == 0L) queue <- c(queue, target)
    }
  }
  if (length(unique(visited)) != length(node_ids)) stop("`edges` must form an acyclic graph.", call. = FALSE)
  layout <- data.frame(ID = node_ids, Level = unname(level[node_ids]), stringsAsFactors = FALSE)
  layout$X <- ave(layout$Level, layout$Level, FUN = function(x) seq_along(x) - (length(x) + 1) / 2)
  layout$Y <- -layout$Level
  node_data <- merge(nodes, layout, by.x = id, by.y = "ID", sort = FALSE)
  node_data <- node_data[match(node_ids, node_data[[id]]), , drop = FALSE]
  node_data$.score <- -log10(node_data[[pvalue]])
  node_data$.label <- paste(node_data[[id]], node_data[[term]],
                            format(node_data[[pvalue]], digits = 2),
                            node_data[[gene_ratio]], sep = "\n")
  edge_data <- merge(edges, layout, by.x = parent, by.y = "ID", sort = FALSE)
  names(edge_data)[names(edge_data) == "X"] <- "XParent"
  names(edge_data)[names(edge_data) == "Y"] <- "YParent"
  edge_data <- merge(edge_data, layout, by.x = child, by.y = "ID", sort = FALSE)
  names(edge_data)[names(edge_data) == "X"] <- "XChild"
  names(edge_data)[names(edge_data) == "Y"] <- "YChild"
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_data, ggplot2::aes(
      x = .data$XParent, y = .data$YParent, xend = .data$XChild, yend = .data$YChild
    ), color = "#777777", linewidth = 0.5,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
    ggplot2::geom_label(data = node_data, ggplot2::aes(
      x = .data$X, y = .data$Y, label = .data$.label, fill = .data$.score
    ), family = z$style$global$font_family,
    size = z$style$text$data_label$size / 3.2, label.size = 0.25) +
    ggplot2::scale_fill_gradient(low = "#FFF5B1", high = "#C51B1D",
                                 name = "-log10(Pvalue)",
                                 breaks = .default_numeric_breaks(node_data$.score, 3)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.28))) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = title) +
    z$style$ggplot_theme +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(
                     family = z$style$global$font_family,
                     size = z$style$text$title$size,
                     face = .text_face(z$style$text$title), margin = ggplot2::margin(b = 14)
                   ))
  .sanshu_finish(p, z)
}

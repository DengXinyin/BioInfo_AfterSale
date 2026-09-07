# Set and network plots consolidated from venn_upset.r, plot_network.r and
# zi_pi.R in micro_v2:v5.42.

micro_overlap_plot <- function(sets, colors = NULL, show_percentage = FALSE,
                               font_family = "sans") {
  if (!is.list(sets) || length(sets) < 2) stop("`sets` must be a named list with at least two sets.")
  if (is.null(names(sets))) names(sets) <- paste0("Set", seq_along(sets))
  sets <- lapply(sets, unique)
  if (length(sets) <= 4 && requireNamespace("ggvenn", quietly = TRUE)) {
    if (is.null(colors)) colors <- .micro_palette(length(sets))
    return(ggvenn::ggvenn(sets, fill_color = colors, stroke_size = 0.5,
                          set_name_size = 4, text_size = 4,
                          show_percentage = show_percentage) +
             ggplot2::theme(text = ggplot2::element_text(family = font_family)))
  }
  stop("Two to four sets and package `ggvenn` are currently required.")
}

micro_network_plot <- function(edges, nodes = NULL, source = "Source", target = "Target",
                               weight = "Weight", node = "Node", node_group = NULL,
                               top_edges = 50, seed = 123, colors = NULL) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package `igraph` is required.")
  if (!all(c(source, target) %in% names(edges))) stop("Source and target columns are required.")
  e <- edges
  if (!weight %in% names(e)) e[[weight]] <- 1
  e <- e[order(abs(e[[weight]]), decreasing = TRUE), , drop = FALSE]
  e <- utils::head(e, top_edges)
  names(e)[match(c(source, target), names(e))] <- c("from", "to")
  vertices <- unique(c(as.character(e$from), as.character(e$to)))
  if (is.null(nodes)) nodes <- data.frame(Node = vertices, stringsAsFactors = FALSE)
  names(nodes)[match(node, names(nodes))] <- "name"
  nodes <- nodes[nodes$name %in% vertices, , drop = FALSE]
  graph <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = nodes)
  set.seed(seed)
  layout <- igraph::layout_with_fr(graph)
  if (!is.null(node_group) && node_group %in% names(nodes)) {
    fac <- factor(igraph::vertex_attr(graph, node_group))
    if (is.null(colors)) colors <- .micro_palette(nlevels(fac))
    vertex_colors <- colors[fac]
  } else vertex_colors <- "#4DAF4A"
  edge_weights <- igraph::edge_attr(graph, weight)
  edge_colors <- ifelse(edge_weights >= 0, "#E69F00", "#56B4E9")
  draw <- function() {
    plot(graph, layout = layout, vertex.color = vertex_colors, vertex.frame.color = NA,
         vertex.size = 10, vertex.label.cex = 0.75,
         edge.color = edge_colors, edge.width = 1 + 3 * abs(edge_weights) /
           max(abs(edge_weights)))
  }
  structure(list(plot = draw, graph = graph, layout = layout), class = "micro_network")
}

micro_zipi_plot <- function(edges, source = "Source", target = "Target", weight = "Weight",
                            node_group = NULL, colors = NULL, font_family = "sans") {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Package `igraph` is required.")
  e <- edges
  names(e)[match(c(source, target), names(e))] <- c("from", "to")
  if (!weight %in% names(e)) e[[weight]] <- 1
  g <- igraph::graph_from_data_frame(e, directed = FALSE)
  membership <- igraph::membership(igraph::cluster_fast_greedy(
    g, weights = abs(igraph::edge_attr(g, weight))))
  a <- as.matrix(igraph::as_adjacency_matrix(g, attr = weight, sparse = FALSE))
  degree <- rowSums(abs(a))
  zi <- numeric(nrow(a)); pi <- numeric(nrow(a))
  for (m in unique(membership)) {
    idx <- membership == m
    within <- rowSums(abs(a[, idx, drop = FALSE]))
    s <- stats::sd(within[idx])
    zi[idx] <- if (is.finite(s) && s > 0) (within[idx] - mean(within[idx])) / s else 0
  }
  for (i in seq_len(nrow(a))) {
    by_module <- vapply(unique(membership), function(m) sum(abs(a[i, membership == m])), numeric(1))
    pi[i] <- if (degree[i] > 0) 1 - sum((by_module / degree[i])^2) else 0
  }
  dat <- data.frame(Node = rownames(a), Zi = zi, Pi = pi, Module = factor(membership))
  dat$Role <- ifelse(dat$Zi > 2.5 & dat$Pi > 0.62, "Network hub",
                     ifelse(dat$Zi > 2.5, "Module hub",
                            ifelse(dat$Pi > 0.62, "Connector", "Peripheral")))
  if (is.null(colors)) colors <- .micro_palette(nlevels(dat$Module))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Pi, y = Zi, colour = Module)) +
    ggplot2::geom_vline(xintercept = 0.62, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = 2.5, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 2.5) + ggplot2::scale_colour_manual(values = colors) +
    ggplot2::labs(x = "Among-module connectivity (Pi)",
                  y = "Within-module connectivity (Zi)", title = "Network node roles") +
    .micro_theme(font_family)
  structure(list(plot = p, nodes = dat, graph = g), class = "micro_zipi")
}

# Binning diagnostics and gene-set overlap plots.

metage_gc_coverage <- function(data, gc = "GC", coverage = "Coverage", bin = "Bin",
                                colors = NULL, font_family = metage_font_family()) {
  if (!all(c(gc, coverage, bin) %in% names(data))) stop("GC/coverage/bin columns are missing.")
  data[[gc]] <- suppressWarnings(as.numeric(data[[gc]]))
  data[[coverage]] <- suppressWarnings(as.numeric(data[[coverage]]))
  dat <- data[stats::complete.cases(data[c(gc, coverage, bin)]), , drop = FALSE]
  dat$.Bin <- factor(dat[[bin]], levels = unique(dat[[bin]]))
  if (is.null(colors)) colors <- .metage_palette(nlevels(dat$.Bin))
  ggplot2::ggplot(dat, ggplot2::aes(.data[[gc]], .data[[coverage]], colour = .Bin)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.5) + ggplot2::scale_colour_manual(values = colors) +
    ggplot2::labs(x = "GC content", y = "Coverage", title = "Bin GC-coverage") +
    .metage_theme(font_family)
}

metage_bin_heatmap <- function(abundance, ...) {
  metage_abundance_heatmap(abundance, title = "Bin abundance", ...)
}

metage_overlap_plot <- function(sets, colors = NULL, show_percentage = FALSE,
                                 font_family = metage_font_family()) {
  if (!is.list(sets) || length(sets) < 2) stop("`sets` must be a list of at least two sets.")
  if (is.null(names(sets))) names(sets) <- paste0("Set", seq_along(sets))
  sets <- lapply(sets, unique)
  if (length(sets) <= 4 && requireNamespace("ggvenn", quietly = TRUE)) {
    if (is.null(colors)) colors <- .metage_palette(length(sets))
    return(ggvenn::ggvenn(sets, fill_color = colors, show_percentage = show_percentage,
                          stroke_size = 0.5, set_name_size = 5, text_size = 4) +
             ggplot2::theme(text = ggplot2::element_text(family = font_family)))
  }
  if (requireNamespace("UpSetR", quietly = TRUE)) {
    return(function() UpSetR::upset(UpSetR::fromList(sets), order.by = "freq"))
  }
  stop("Install `ggvenn` (2-4 sets) or `UpSetR` (>4 sets).")
}

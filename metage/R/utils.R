# Shared helpers for visualizations reconstructed from metage:v2.87.

`%metage_or%` <- function(x, y) if (is.null(x)) y else x

.metage_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package `", package, "` is required.", call. = FALSE)
  }
}

.metage_matrix <- function(x, name = "abundance") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", name, "` must be a numeric matrix/data frame (features x samples).",
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (!nrow(x) || !ncol(x)) stop("`", name, "` must not be empty.", call. = FALSE)
  if (is.null(rownames(x))) rownames(x) <- paste0("Feature_", seq_len(nrow(x)))
  if (is.null(colnames(x))) colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` contains missing/non-finite values.", call. = FALSE)
  }
  x
}

.metage_group <- function(group, samples) {
  if (is.null(group)) return(stats::setNames(factor(samples), samples))
  if (!is.null(names(group))) {
    absent <- setdiff(samples, names(group))
    if (length(absent)) stop("`group` is missing: ", paste(absent, collapse = ", "))
    group <- group[samples]
  } else if (length(group) != length(samples)) {
    stop("Unnamed `group` must contain one value per sample.", call. = FALSE)
  }
  out <- factor(as.character(group), levels = unique(as.character(group)))
  names(out) <- samples
  out
}

.metage_palette <- function(n) {
  base <- c("#377EB8", "#FF7F00", "#4DAF4A", "#E41A1C", "#984EA3",
            "#A65628", "#F781BF", "#20B2AA", "#CE50CA", "#CBD588")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

metage_font_family <- function(preferred = c("Times New Roman", "Arial")) {
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    available <- unique(systemfonts::system_fonts()$family)
    hit <- preferred[preferred %in% available]
    if (length(hit)) return(hit[1])
  }
  "sans"
}

.metage_theme <- function(font_family = metage_font_family(), base_size = 16) {
  ggplot2::theme_bw(base_family = font_family, base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 20, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 16, colour = "black"),
      legend.title = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 14)
    )
}

.metage_long <- function(x, feature = "Feature", sample = "Sample") {
  out <- data.frame(
    Feature = rep(rownames(x), times = ncol(x)),
    Sample = rep(colnames(x), each = nrow(x)),
    Value = as.vector(x), check.names = FALSE, stringsAsFactors = FALSE
  )
  names(out)[1:2] <- c(feature, sample)
  out
}

metage_save_plot <- function(plot, output_file, width = 8, height = 6, dpi = 300) {
  stopifnot(is.character(output_file), length(output_file) == 1L)
  ext <- tolower(tools::file_ext(output_file))
  if (!ext %in% c("pdf", "png")) stop("Output must end in .pdf or .png.")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(output_file, plot, width = width, height = height, dpi = dpi,
                    device = if (ext == "pdf") grDevices::cairo_pdf else "png")
  } else {
    if (ext == "pdf") grDevices::cairo_pdf(output_file, width, height)
    else grDevices::png(output_file, width = width, height = height, units = "in", res = dpi)
    on.exit(grDevices::dev.off(), add = TRUE)
    if (inherits(plot, c("grob", "gTree", "gtable"))) {
      grid::grid.newpage(); grid::grid.draw(plot)
    } else if (is.function(plot)) plot() else print(plot)
  }
  invisible(normalizePath(output_file, mustWork = TRUE))
}

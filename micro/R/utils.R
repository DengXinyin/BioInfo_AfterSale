# Internal helpers shared by the micro_v2 visualization functions.

`%||%` <- function(x, y) if (is.null(x)) y else x

.micro_palette <- function(n) {
  base <- c(
    "#377EB8", "#FF7F00", "#4DAF4A", "#E41A1C", "#984EA3",
    "#A65628", "#F781BF", "#20B2AA", "#CBD588", "#CE50CA"
  )
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

.micro_matrix <- function(x, name = "abundance") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", name, "` must be a numeric matrix or data frame.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (is.null(rownames(x))) rownames(x) <- paste0("Feature_", seq_len(nrow(x)))
  if (is.null(colnames(x))) colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` contains missing or non-finite values.", call. = FALSE)
  }
  x
}

.micro_group <- function(group, sample_names) {
  if (is.null(group)) group <- stats::setNames(sample_names, sample_names)
  if (length(group) != length(sample_names) && is.null(names(group))) {
    stop("Unnamed `group` must have one value per sample.", call. = FALSE)
  }
  if (!is.null(names(group))) {
    missing <- setdiff(sample_names, names(group))
    if (length(missing)) {
      stop("`group` is missing samples: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    group <- group[sample_names]
  }
  out <- factor(as.character(group), levels = unique(as.character(group)))
  names(out) <- sample_names
  out
}

.micro_long_matrix <- function(x, row_name = "Feature", col_name = "Sample") {
  out <- data.frame(
    Feature = rep(rownames(x), times = ncol(x)),
    Sample = rep(colnames(x), each = nrow(x)),
    Value = as.vector(x),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stats::setNames(out, c(row_name, col_name, "Value"))
}

.micro_theme <- function(font_family = "sans", base_size = 12) {
  ggplot2::theme_bw(base_family = font_family, base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_blank()
    )
}

micro_save_plot <- function(plot, output_file, width = 8, height = 6, dpi = 300) {
  if (!is.character(output_file) || length(output_file) != 1L) {
    stop("`output_file` must be one file path.", call. = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(output_file))
  if (!ext %in% c("png", "pdf")) stop("Output must end in .png or .pdf.", call. = FALSE)

  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(output_file, plot, width = width, height = height, dpi = dpi,
                    device = if (ext == "pdf") grDevices::cairo_pdf else "png")
  } else {
    if (ext == "png") grDevices::png(output_file, width = width, height = height,
                                      units = "in", res = dpi)
    else grDevices::cairo_pdf(output_file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    if (inherits(plot, c("grob", "gTree", "gtable"))) {
      grid::grid.newpage()
      grid::grid.draw(plot)
    } else if (inherits(plot, "recordedplot")) {
      grDevices::replayPlot(plot)
    } else if (is.function(plot)) {
      plot()
    } else {
      print(plot)
    }
  }
  invisible(normalizePath(output_file, mustWork = TRUE))
}

.sanshu_style_output <- function(style, output_file, width, height, dpi) {
  resolved <- .validate_style_output(style, output_file, width, height, dpi)
  if (!is.null(output_file) && !grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
    warning("PDF is the recommended output format for transcriptome tutorials.", call. = FALSE)
  }
  resolved
}

.sanshu_require_columns <- function(data, columns) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(data)
}

.sanshu_numeric <- function(data, columns) {
  bad <- columns[!vapply(data[columns], is.numeric, logical(1))]
  if (length(bad)) {
    stop("Column(s) must be numeric: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (any(!is.finite(as.matrix(data[columns])))) {
    stop("Numeric plotting columns must contain only finite values.", call. = FALSE)
  }
  invisible(data)
}

.sanshu_palette <- function(levels, colors, style) {
  if (is.null(colors)) colors <- rep(style$group_palette, length.out = length(levels))
  if (is.null(names(colors))) names(colors) <- levels
  missing <- setdiff(levels, names(colors))
  if (length(missing)) {
    stop("Colors are missing level(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  colors[levels]
}

.sanshu_finish <- function(plot, resolved) {
  .save_bioinfo_plot(
    plot, resolved$style, resolved$output_file,
    resolved$width, resolved$height, resolved$dpi
  )
  plot
}

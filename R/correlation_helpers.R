# Shared validation and correlation helpers for correlation/RDA plots.
.correlation_matrix <- function(x, y, method = "spearman") {
  x <- as.data.frame(x, check.names = FALSE); y <- as.data.frame(y, check.names = FALSE)
  if (!nrow(x) || !nrow(y)) stop("Both inputs must contain rows.", call. = FALSE)
  if (is.null(rownames(x)) || is.null(rownames(y))) stop("Both inputs must have row names identifying samples.", call. = FALSE)
  common <- intersect(rownames(x), rownames(y)); if (length(common) < 3L) stop("At least three shared samples are required.", call. = FALSE)
  x <- x[common, , drop = FALSE]; y <- y[common, , drop = FALSE]; x[] <- lapply(x, as.numeric); y[] <- lapply(y, as.numeric)
  out <- matrix(NA_real_, nrow = ncol(x), ncol = ncol(y), dimnames = list(colnames(x), colnames(y)))
  for (i in seq_len(ncol(x))) for (j in seq_len(ncol(y))) {
    ok <- is.finite(x[[i]]) & is.finite(y[[j]])
    if (sum(ok) >= 3L && length(unique(x[[i]][ok])) > 1L && length(unique(y[[j]][ok])) > 1L) out[i, j] <- suppressWarnings(stats::cor(x[[i]][ok], y[[j]][ok], method = method))
  }
  out
}
.correlation_long <- function(rho) { d <- as.data.frame(as.table(rho), stringsAsFactors = FALSE); names(d) <- c("row", "column", "value"); d$size <- abs(d$value); d }
.validate_style_output <- function(style, output_file, width, height, dpi) {
  style <- .plot_style_or_default(style)
  list(style = style, output_file = output_file, width = width %||% style$global$figure_width, height = height %||% style$global$figure_height, dpi = dpi %||% style$global$dpi)
}

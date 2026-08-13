#' @export
print.bioinfo_enrichment <- function(x, ...) {
  cat("BioInfoAfterSale enrichment result\n")
  cat("  Terms:", nrow(x$table), "\n")
  cat("  Components: result, table, plot, call\n")
  invisible(x)
}

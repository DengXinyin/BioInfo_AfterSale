.plot_null_default <- function(left, right) {
  if (is.null(left)) right else left
}

`%||%` <- .plot_null_default

.assert_choice <- function(value, choices, argument) {
  if (length(value) != 1L || is.na(value) || !value %in% choices) {
    stop(
      sprintf("`%s` must be one of: %s.", argument, paste(choices, collapse = ", ")),
      call. = FALSE
    )
  }
  value
}

.assert_probability <- function(value, argument) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      value < 0 || value > 1) {
    stop(sprintf("`%s` must be a number between 0 and 1.", argument), call. = FALSE)
  }
  value
}

.assert_positive_number <- function(value, argument) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) || value <= 0) {
    stop(sprintf("`%s` must be a positive number.", argument), call. = FALSE)
  }
  value
}

.assert_numeric_vector <- function(
    value, argument, length_required = NULL, positive = FALSE) {
  if (!is.numeric(value) || !length(value) || anyNA(value) ||
      any(!is.finite(value))) {
    stop(sprintf("`%s` must be a finite numeric vector.", argument), call. = FALSE)
  }
  if (!is.null(length_required) && length(value) != length_required) {
    stop(
      sprintf("`%s` must contain exactly %d values.", argument, length_required),
      call. = FALSE
    )
  }
  if (positive && any(value <= 0)) {
    stop(sprintf("`%s` values must be positive.", argument), call. = FALSE)
  }
  value
}

.assert_flag <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", argument), call. = FALSE)
  }
  value
}

.default_size_breaks <- function(values, n = 3L) {
  values <- sort(unique(values[is.finite(values)]))
  if (length(values) <= n) return(values)
  positions <- unique(floor(seq(1, length(values), length.out = n)))
  values[positions]
}

.default_numeric_breaks <- function(values, n = 3L) {
  values <- values[is.finite(values)]
  if (!length(values)) return(numeric())
  if (diff(range(values)) == 0) return(values[1])
  unique(signif(seq(min(values), max(values), length.out = n), 3))
}

.escape_html <- function(values) {
  values <- gsub("&", "&amp;", values, fixed = TRUE)
  values <- gsub("<", "&lt;", values, fixed = TRUE)
  gsub(">", "&gt;", values, fixed = TRUE)
}

.clean_gene_ids <- function(gene, argument = "gene") {
  if (is.factor(gene)) gene <- as.character(gene)
  if (!is.atomic(gene)) {
    stop(sprintf("`%s` must be an atomic vector of gene IDs.", argument), call. = FALSE)
  }
  gene <- trimws(as.character(gene))
  gene <- unique(gene[!is.na(gene) & nzchar(gene)])
  if (!length(gene)) {
    stop(sprintf("`%s` contains no usable gene IDs.", argument), call. = FALSE)
  }
  gene
}

.resolve_orgdb <- function(org_db = NULL, species = NULL) {
  if (!is.null(org_db)) return(org_db)
  if (is.null(species)) {
    stop("GO analysis requires either `org_db` or `species`.", call. = FALSE)
  }

  species_key <- tolower(trimws(species))
  aliases <- c(
    human = "org.Hs.eg.db", homo_sapiens = "org.Hs.eg.db", hsa = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db", mus_musculus = "org.Mm.eg.db", mmu = "org.Mm.eg.db",
    rat = "org.Rn.eg.db", rattus_norvegicus = "org.Rn.eg.db", rno = "org.Rn.eg.db"
  )
  package <- unname(aliases[[species_key]])
  if (is.null(package)) {
    stop(
      "Unsupported `species`. Use human/hsa, mouse/mmu, rat/rno, or provide `org_db`.",
      call. = FALSE
    )
  }
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this GO analysis.", package), call. = FALSE)
  }
  get(package, envir = asNamespace(package))
}

.normalise_result_table <- function(result) {
  if (methods::is(result, "enrichResult") || methods::is(result, "gseaResult")) {
    return(as.data.frame(result))
  }
  if (is.data.frame(result)) return(result)
  stop("`result` must be an enrichResult, gseaResult, or data.frame.", call. = FALSE)
}

.filter_result_object <- function(result, filter_by, cutoff) {
  if (!methods::is(result, "enrichResult") && !methods::is(result, "gseaResult")) {
    stop("Filtering for enrichplot requires an enrichResult or gseaResult object.", call. = FALSE)
  }
  table <- as.data.frame(result)
  if (!filter_by %in% names(table)) {
    stop(sprintf("Column '%s' is not present in the enrichment result.", filter_by), call. = FALSE)
  }
  values <- suppressWarnings(as.numeric(table[[filter_by]]))
  zero_count <- sum(!is.na(values) & values == 0)
  if (zero_count) {
    warning(
      sprintf(
        "%d row(s) with %s = 0 were excluded; valid enrichment filtering requires 0 < p < cutoff.",
        zero_count, filter_by
      ),
      call. = FALSE
    )
  }
  keep <- !is.na(values) & values > 0 & values < cutoff
  result@result <- result@result[keep, , drop = FALSE]
  result
}

.filter_positive_pvalues <- function(values, filter_by, cutoff) {
  values <- suppressWarnings(as.numeric(values))
  zero_count <- sum(!is.na(values) & values == 0)
  if (zero_count) {
    warning(
      sprintf(
        "%d row(s) with %s = 0 were excluded; valid enrichment filtering requires 0 < p < cutoff.",
        zero_count, filter_by
      ),
      call. = FALSE
    )
  }
  !is.na(values) & values > 0 & values < cutoff
}
# Plotting functions in this package prefer Times New Roman for all text.

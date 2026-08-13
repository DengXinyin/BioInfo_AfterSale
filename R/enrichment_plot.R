#' Run enrichment analysis and draw the result in one call
#'
#' This is the main convenience interface. Supply `gene` to run a new GO or
#' KEGG analysis, or supply an existing `enrichResult`/`gseaResult` through
#' `result` to skip calculation and only redraw the figure.
#'
#' @param gene Optional character vector of gene IDs.
#' @param result Optional existing enrichment result object.
#' @param analysis Either `"GO"` or `"KEGG"`; used when `gene` is supplied.
#' @param run_args Named list passed to [run_enrichment()].
#' @param plot_args Named list passed to [plot_enrichment()].
#'
#' @return An object of class `bioinfo_enrichment` containing `result`,
#'   `table`, `plot`, and `call`.
#' @export
enrichment_plot <- function(
    gene = NULL,
    result = NULL,
    analysis = c("GO", "KEGG"),
    run_args = list(),
    plot_args = list()) {

  has_gene <- !is.null(gene)
  has_result <- !is.null(result)
  if (has_gene == has_result) {
    stop("Supply exactly one of `gene` or `result`.", call. = FALSE)
  }
  if (!is.list(run_args)) stop("`run_args` must be a named list.", call. = FALSE)
  if (!is.list(plot_args)) stop("`plot_args` must be a named list.", call. = FALSE)

  if (has_gene) {
    analysis <- match.arg(analysis)
    protected <- intersect(names(run_args), c("gene", "analysis"))
    if (length(protected)) {
      stop("Do not include `gene` or `analysis` inside `run_args`.", call. = FALSE)
    }
    result <- do.call(
      run_enrichment,
      c(list(gene = gene, analysis = analysis), run_args)
    )
  }

  if ("result" %in% names(plot_args)) {
    stop("Do not include `result` inside `plot_args`.", call. = FALSE)
  }
  plot <- do.call(plot_enrichment, c(list(result = result), plot_args))
  output <- list(
    result = result,
    table = .normalise_result_table(result),
    plot = plot,
    call = match.call()
  )
  class(output) <- c("bioinfo_enrichment", "list")
  output
}

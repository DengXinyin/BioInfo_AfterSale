#' Run GO or KEGG enrichment analysis
#'
#' Runs over-representation analysis from a vector of gene IDs. GO analysis
#' uses [clusterProfiler::enrichGO()]. KEGG analysis can use either a local
#' `TERM2GENE` table via [clusterProfiler::enricher()] or KEGG through
#' [clusterProfiler::enrichKEGG()]. Local annotation is recommended in a
#' reproducible or offline Docker environment.
#'
#' @param gene Character vector of input gene IDs.
#' @param analysis Either `"GO"` or `"KEGG"`.
#' @param universe Optional character vector of background gene IDs.
#' @param species Species alias for GO (`human`, `mouse`, or `rat`). For online
#'   KEGG analysis this is a KEGG organism code such as `hsa` or `mmu`.
#' @param key_type Gene ID type. The default is `"ENTREZID"` for GO and
#'   `"kegg"` for online KEGG. It can be set to `"SYMBOL"` or `"ENSEMBL"`
#'   for GO when supported by the selected OrgDb.
#' @param org_db Optional OrgDb object. When supplied it takes precedence over
#'   `species` for GO analysis.
#' @param ont GO ontology: `"ALL"`, `"BP"`, `"CC"`, or `"MF"`.
#' @param term2gene Optional two-column pathway-to-gene data frame for local
#'   KEGG analysis.
#' @param term2name Optional two-column pathway-to-name data frame.
#' @param pvalue_cutoff P-value cutoff passed to clusterProfiler.
#' @param p_adjust_method Multiple-testing correction method.
#' @param qvalue_cutoff Q-value cutoff passed to clusterProfiler.
#' @param min_gs_size Minimum gene-set size.
#' @param max_gs_size Maximum gene-set size.
#' @param readable Convert GO result gene IDs to symbols when possible.
#' @param ... Additional arguments passed to the underlying clusterProfiler
#'   function.
#'
#' @return An `enrichResult` object.
#' @export
GO_KEGG_analyse <- function(
    gene,
    analysis = c("GO", "KEGG"),
    universe = NULL,
    species = NULL,
    key_type = NULL,
    org_db = NULL,
    ont = "ALL",
    term2gene = NULL,
    term2name = NULL,
    pvalue_cutoff = 0.05,
    p_adjust_method = "BH",
    qvalue_cutoff = 0.2,
    min_gs_size = 10,
    max_gs_size = 500,
    readable = TRUE,
    ...) {

  analysis <- match.arg(analysis)
  gene <- .clean_gene_ids(gene)
  if (!is.null(universe)) universe <- .clean_gene_ids(universe, "universe")
  .assert_probability(pvalue_cutoff, "pvalue_cutoff")
  .assert_probability(qvalue_cutoff, "qvalue_cutoff")

  if (analysis == "GO") {
    if (is.null(key_type)) key_type <- "ENTREZID"
    ont <- .assert_choice(toupper(ont), c("ALL", "BP", "CC", "MF"), "ont")
    org_db <- .resolve_orgdb(org_db, species)
    return(clusterProfiler::enrichGO(
      gene = gene, universe = universe, OrgDb = org_db, keyType = key_type,
      ont = ont, pAdjustMethod = p_adjust_method, pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff, minGSSize = min_gs_size,
      maxGSSize = max_gs_size, readable = readable, ...
    ))
  }

  if (!is.null(term2gene)) {
    if (!is.data.frame(term2gene) || ncol(term2gene) < 2L) {
      stop("`term2gene` must be a data frame with at least two columns.", call. = FALSE)
    }
    if (!is.null(term2name) && (!is.data.frame(term2name) || ncol(term2name) < 2L)) {
      stop("`term2name` must be a data frame with at least two columns.", call. = FALSE)
    }
    return(clusterProfiler::enricher(
      gene = gene, universe = universe,
      TERM2GENE = term2gene[, 1:2, drop = FALSE],
      TERM2NAME = if (is.null(term2name)) NULL else term2name[, 1:2, drop = FALSE],
      pAdjustMethod = p_adjust_method, pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff, minGSSize = min_gs_size,
      maxGSSize = max_gs_size, ...
    ))
  }

  if (is.null(species)) {
    stop("KEGG analysis requires `species` or a local `term2gene` table.", call. = FALSE)
  }
  if (is.null(key_type)) key_type <- "kegg"
  clusterProfiler::enrichKEGG(
    gene = gene, universe = universe, organism = species, keyType = key_type,
    pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_adjust_method,
    qvalueCutoff = qvalue_cutoff, minGSSize = min_gs_size,
    maxGSSize = max_gs_size, ...
  )
}
# All downstream GO/KEGG plot text prefers Times New Roman.

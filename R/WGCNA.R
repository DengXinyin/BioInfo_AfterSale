# Dependency policy
# -----------------
# Do not install packages from inside analysis functions. Package installation
# belongs to environment/image construction. Each public function below lists
# the R packages it uses and calls non-base dependencies with `package::fun`.
# The only temporary attachment is WGCNA inside WGCNA_build_modules(), required
# by WGCNA 1.73's internal Pearson-correlation lookup; it is detached on exit.

#' Prepare an expression matrix for WGCNA
#'
#' Validates a gene-by-sample expression matrix, optionally transforms it,
#' applies the zero-expression rule used by the reference workflow, transposes
#' it to the sample-by-gene layout required by WGCNA, and runs
#' [WGCNA::goodSamplesGenes()]. Set `keep_all_genes = TRUE` when a customer
#' explicitly requires every supplied gene to enter network construction.
#'
#' @param expression Numeric matrix or data frame. Genes are rows and samples
#'   are columns unless `samples_in_rows = TRUE`.
#' @param gene_id_col Optional column name or index containing gene IDs. It is
#'   removed from the numeric matrix and used as row names.
#' @param samples_in_rows Whether `expression` is already sample-by-gene.
#' @param transform Expression transformation: `"none"` or `"log2p1"`.
#' @param zero_fraction_max Maximum permitted fraction of zero values per gene.
#'   Ignored when `keep_all_genes = TRUE`; use `NULL` to disable this filter.
#' @param keep_all_genes Disable the zero-fraction filter.
#' @param remove_bad Whether samples or genes rejected by
#'   [WGCNA::goodSamplesGenes()] should be removed. If `FALSE`, the function
#'   stops and reports them.
#' @param verbose Verbosity passed to [WGCNA::goodSamplesGenes()].
#'
#' @return A `WGCNA_prepared` list containing `datExpr`, the gene-by-sample
#'   matrix, QC statistics, and rejected IDs.
#' @export
WGCNA_prepare_expression <- function(
    expression,
    gene_id_col = NULL,
    samples_in_rows = FALSE,
    transform = c("none", "log2p1"),
    zero_fraction_max = 0.3,
    keep_all_genes = FALSE,
    remove_bad = TRUE,
    verbose = 0) {

  # Required R packages: WGCNA
  # Base/recommended packages: stats
  WGCNA_require_namespace()
  transform <- match.arg(transform)
  WGCNA_assert_flag(samples_in_rows, "samples_in_rows")
  WGCNA_assert_flag(keep_all_genes, "keep_all_genes")
  WGCNA_assert_flag(remove_bad, "remove_bad")

  x <- as.data.frame(expression, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(gene_id_col)) {
    if (samples_in_rows) {
      stop("`gene_id_col` cannot be used when `samples_in_rows = TRUE`.", call. = FALSE)
    }
    if (is.character(gene_id_col)) {
      if (length(gene_id_col) != 1L || !gene_id_col %in% names(x)) {
        stop("`gene_id_col` is not present in `expression`.", call. = FALSE)
      }
      gene_index <- match(gene_id_col, names(x))
    } else if (is.numeric(gene_id_col) && length(gene_id_col) == 1L &&
               is.finite(gene_id_col) && gene_id_col == as.integer(gene_id_col) &&
               gene_id_col >= 1L && gene_id_col <= ncol(x)) {
      gene_index <- as.integer(gene_id_col)
    } else {
      stop("`gene_id_col` must be one valid column name or index.", call. = FALSE)
    }
    gene_ids <- trimws(as.character(x[[gene_index]]))
    x <- x[, -gene_index, drop = FALSE]
    rownames(x) <- gene_ids
  }

  if (!nrow(x) || !ncol(x)) stop("`expression` must not be empty.", call. = FALSE)
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop("Expression gene/sample row names must be non-empty.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) stop("Expression row names must be unique.", call. = FALSE)
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
    stop("Expression column names must be non-empty.", call. = FALSE)
  }
  if (anyDuplicated(colnames(x))) stop("Expression column names must be unique.", call. = FALSE)

  original_na <- is.na(x)
  numeric_x <- suppressWarnings(as.data.frame(
    lapply(x, as.numeric), check.names = FALSE, stringsAsFactors = FALSE
  ))
  names(numeric_x) <- names(x)
  rownames(numeric_x) <- rownames(x)
  introduced_na <- is.na(numeric_x) & !original_na
  if (any(introduced_na)) {
    stop("`expression` contains non-numeric values.", call. = FALSE)
  }
  matrix_x <- as.matrix(numeric_x)
  if (anyNA(matrix_x) || any(!is.finite(matrix_x))) {
    stop("`expression` contains missing or non-finite values.", call. = FALSE)
  }
  if (transform == "log2p1") {
    if (any(matrix_x < 0)) {
      stop("`log2p1` requires non-negative expression values.", call. = FALSE)
    }
    matrix_x <- log2(matrix_x + 1)
  }

  if (samples_in_rows) {
    datExpr <- as.data.frame(matrix_x, check.names = FALSE)
    gene_sample <- t(matrix_x)
  } else {
    gene_sample <- matrix_x
    datExpr <- as.data.frame(t(matrix_x), check.names = FALSE)
  }

  initial_genes <- ncol(datExpr)
  initial_samples <- nrow(datExpr)
  zero_fraction <- colMeans(datExpr == 0)
  zero_removed <- character()
  if (!keep_all_genes && !is.null(zero_fraction_max)) {
    WGCNA_assert_probability(zero_fraction_max, "zero_fraction_max")
    keep <- zero_fraction <= zero_fraction_max
    zero_removed <- colnames(datExpr)[!keep]
    datExpr <- datExpr[, keep, drop = FALSE]
  }
  if (!ncol(datExpr)) stop("No genes remain after zero-expression filtering.", call. = FALSE)

  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = verbose)
  bad_genes <- names(datExpr)[!gsg$goodGenes]
  bad_samples <- rownames(datExpr)[!gsg$goodSamples]
  if (!gsg$allOK && !remove_bad) {
    stop(
      "goodSamplesGenes rejected genes [", paste(bad_genes, collapse = ", "),
      "] and samples [", paste(bad_samples, collapse = ", "), "].",
      call. = FALSE
    )
  }
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  if (nrow(datExpr) < 4L || ncol(datExpr) < 2L) {
    stop("WGCNA requires at least 4 usable samples and 2 usable genes.", call. = FALSE)
  }

  gene_sample <- t(as.matrix(datExpr))
  result <- list(
    datExpr = datExpr,
    expression = gene_sample,
    qc = data.frame(
      metric = c(
        "input_genes", "input_samples", "output_genes", "output_samples",
        "zero_filter_removed_genes", "goodSamplesGenes_removed_genes",
        "goodSamplesGenes_removed_samples"
      ),
      value = c(
        initial_genes, initial_samples, ncol(datExpr), nrow(datExpr),
        length(zero_removed), length(bad_genes), length(bad_samples)
      ),
      stringsAsFactors = FALSE
    ),
    zero_fraction = zero_fraction,
    removed = list(
      zero_filter_genes = zero_removed,
      bad_genes = bad_genes,
      bad_samples = bad_samples
    ),
    parameters = list(
      transform = transform,
      zero_fraction_max = zero_fraction_max,
      keep_all_genes = keep_all_genes,
      remove_bad = remove_bad
    ),
    sample_tree = stats::hclust(stats::dist(datExpr), method = "average")
  )
  class(result) <- c("WGCNA_prepared", "list")
  result
}

#' Select a WGCNA soft-thresholding power
#'
#' Evaluates scale-free topology over candidate powers. `strategy =
#' "scale_free"` selects the first candidate reaching `fit_cutoff` with a
#' negative slope and falls back to the sample-count rule from the reference
#' image. `strategy = "image_rule"` reproduces that rule directly.
#'
#' @param prepared A `WGCNA_prepared` object or numeric sample-by-gene matrix.
#' @param powers Candidate powers.
#' @param network_type WGCNA network type.
#' @param strategy Selection strategy: `"scale_free"`, `"image_rule"`, or
#'   `"manual"`.
#' @param fit_cutoff Required signed scale-free topology fit.
#' @param manual_power Power used by `strategy = "manual"`.
#' @param verbose Verbosity passed to [WGCNA::pickSoftThreshold()].
#'
#' @return A `WGCNA_power` list containing fit indices, selected power, and
#'   selection reason.
#' @export
WGCNA_select_power <- function(
    prepared,
    powers = c(1:10, seq(12, 20, by = 2)),
    network_type = c("unsigned", "signed", "signed hybrid"),
    strategy = c("scale_free", "image_rule", "manual"),
    fit_cutoff = 0.9,
    manual_power = NULL,
    verbose = 0) {

  # Required R packages: WGCNA
  # Base packages: stats
  WGCNA_require_namespace()
  datExpr <- WGCNA_get_datExpr(prepared)
  network_type <- match.arg(network_type)
  strategy <- match.arg(strategy)
  WGCNA_assert_probability(fit_cutoff, "fit_cutoff")
  if (!is.numeric(powers) || !length(powers) || anyNA(powers) ||
      any(!is.finite(powers)) || any(powers <= 0)) {
    stop("`powers` must be a non-empty vector of positive numbers.", call. = FALSE)
  }

  sft <- WGCNA::pickSoftThreshold(
    datExpr, powerVector = powers, networkType = network_type, verbose = verbose
  )
  fit <- as.data.frame(sft$fitIndices)
  signed_fit <- -sign(fit[[3]]) * fit[[2]]
  image_power <- WGCNA_image_power(nrow(datExpr), network_type)

  if (strategy == "manual") {
    if (is.null(manual_power) || length(manual_power) != 1L ||
        !is.numeric(manual_power) || is.na(manual_power) || manual_power <= 0) {
      stop("`manual_power` must be one positive number.", call. = FALSE)
    }
    selected <- manual_power
    reason <- "manual"
  } else if (strategy == "image_rule") {
    selected <- image_power
    reason <- "reference image sample-count rule"
  } else {
    eligible <- which(is.finite(signed_fit) & signed_fit >= fit_cutoff & fit[[3]] < 0)
    if (length(eligible)) {
      selected <- fit[[1]][eligible[[1]]]
      reason <- sprintf("first candidate with signed R^2 >= %.3f and negative slope", fit_cutoff)
    } else {
      selected <- image_power
      reason <- "no candidate reached fit cutoff; reference image fallback"
    }
  }

  result <- list(
    fit_indices = fit,
    signed_fit = signed_fit,
    power = as.numeric(selected),
    reason = reason,
    network_type = network_type,
    strategy = strategy,
    fit_cutoff = fit_cutoff,
    image_rule_power = image_power
  )
  class(result) <- c("WGCNA_power", "list")
  result
}

#' Build WGCNA co-expression modules
#'
#' Wraps [WGCNA::blockwiseModules()] with the defaults used by the v1.6.7
#' reference image while keeping all important choices configurable.
#'
#' @param prepared A `WGCNA_prepared` object or sample-by-gene matrix.
#' @param power Numeric power or a `WGCNA_power` object.
#' @param network_type,TOM_type Network and TOM type. When `network_type` is
#'   `NULL`, it is inherited from a `WGCNA_power` object or defaults to
#'   `"unsigned"`.
#' @param deep_split,min_module_size,merge_cut_height,reassign_threshold,
#'   pam_respects_dendro,max_block_size,save_TOMs Parameters passed to
#'   [WGCNA::blockwiseModules()].
#' @param threads Number of WGCNA threads.
#' @param seed Random seed.
#' @param verbose WGCNA verbosity.
#'
#' @return A `WGCNA_modules` object containing the WGCNA network, module table,
#'   module colors, and raw eigengenes.
#' @export
WGCNA_build_modules <- function(
    prepared,
    power,
    network_type = NULL,
    TOM_type = NULL,
    deep_split = 2,
    min_module_size = 30,
    merge_cut_height = 0.25,
    reassign_threshold = 0,
    pam_respects_dendro = FALSE,
    max_block_size = 1000,
    save_TOMs = FALSE,
    threads = 1,
    seed = 123,
    verbose = 0) {

  # Required R packages: WGCNA
  # WGCNA imports used internally: dynamicTreeCut, fastcluster
  # Base packages: stats
  WGCNA_require_namespace()
  datExpr <- WGCNA_get_datExpr(prepared)
  power_network_type <- NULL
  if (inherits(power, "WGCNA_power")) {
    power_network_type <- power$network_type
    power <- power$power
  }
  if (is.null(network_type)) network_type <- power_network_type
  if (is.null(network_type)) network_type <- "unsigned"
  network_type <- match.arg(network_type, c("unsigned", "signed", "signed hybrid"))
  WGCNA_assert_positive(power, "power")
  WGCNA_assert_positive(min_module_size, "min_module_size")
  WGCNA_assert_positive(max_block_size, "max_block_size")
  WGCNA_assert_positive(threads, "threads")
  WGCNA_assert_probability(merge_cut_height, "merge_cut_height")
  WGCNA_assert_flag(pam_respects_dendro, "pam_respects_dendro")
  WGCNA_assert_flag(save_TOMs, "save_TOMs")
  if (is.null(TOM_type)) {
    TOM_type <- if (identical(network_type, "signed hybrid")) "signed" else network_type
  }

  # WGCNA 1.73 resolves its Pearson correlation function by name inside
  # blockwiseModules(). When the namespace is loaded with `::` but not
  # attached, that lookup can incorrectly fall through to stats::cor(), which
  # does not accept WGCNA's weight arguments. Attach only for this call and
  # restore the caller's search path afterwards.
  attached_here <- FALSE
  if (!"package:WGCNA" %in% search()) {
    suppressPackageStartupMessages(
      base::library("WGCNA", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
    attached_here <- TRUE
  }
  if (attached_here) {
    on.exit(base::detach("package:WGCNA", unload = FALSE, character.only = TRUE), add = TRUE)
  }

  set.seed(seed)
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    networkType = network_type,
    TOMType = TOM_type,
    deepSplit = deep_split,
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    reassignThreshold = reassign_threshold,
    pamRespectsDendro = pam_respects_dendro,
    numericLabels = FALSE,
    maxBlockSize = max_block_size,
    saveTOMs = save_TOMs,
    nThreads = as.integer(threads),
    verbose = verbose
  )
  colors <- net$colors
  names(colors) <- colnames(datExpr)
  MEs <- WGCNA::orderMEs(net$MEs)
  result <- list(
    net = net,
    datExpr = datExpr,
    module_colors = colors,
    module_gene = data.frame(
      GeneID = colnames(datExpr), module = unname(colors), stringsAsFactors = FALSE
    ),
    MEs_raw = MEs,
    module_summary = data.frame(
      module = names(sort(table(colors), decreasing = TRUE)),
      gene_count = as.integer(sort(table(colors), decreasing = TRUE)),
      stringsAsFactors = FALSE
    ),
    parameters = list(
      power = power, network_type = network_type, TOM_type = TOM_type,
      deep_split = deep_split, min_module_size = min_module_size,
      merge_cut_height = merge_cut_height,
      reassign_threshold = reassign_threshold,
      pam_respects_dendro = pam_respects_dendro,
      max_block_size = max_block_size, threads = as.integer(threads), seed = seed
    )
  )
  class(result) <- c("WGCNA_modules", "list")
  result
}

#' Summarize module eigengenes across samples
#'
#' Produces the statistically meaningful no-trait representation of
#' module-sample relationships. It does not create one-sample dummy traits or
#' report pseudo-correlation P values.
#'
#' @param modules A `WGCNA_modules` object.
#' @param sample_info Optional sample metadata with row names or a `Sample_ID`
#'   column matching expression samples.
#' @param time_col Optional metadata column used for descriptive mean, SD, and
#'   sample count summaries.
#' @param orient Whether each eigengene should be oriented to correlate
#'   positively with its module mean expression.
#'
#' @return A list containing raw/oriented eigengenes, orientation metadata,
#'   long-form sample values, and optional time summaries.
#' @export
WGCNA_module_sample <- function(
    modules,
    sample_info = NULL,
    time_col = NULL,
    orient = TRUE) {

  # Required R packages: none beyond the WGCNA module object
  # Base/recommended packages: stats
  if (!inherits(modules, "WGCNA_modules")) {
    stop("`modules` must be returned by WGCNA_build_modules().", call. = FALSE)
  }
  WGCNA_assert_flag(orient, "orient")
  raw <- modules$MEs_raw
  oriented <- raw
  module_names <- substring(colnames(raw), 3)
  orientation <- data.frame(
    module = module_names,
    raw_ME_vs_module_mean_cor = NA_real_,
    sign_multiplier = 1L,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(module_names)) {
    genes <- modules$module_colors == module_names[[i]]
    module_mean <- rowMeans(modules$datExpr[, genes, drop = FALSE])
    orientation$raw_ME_vs_module_mean_cor[[i]] <- stats::cor(raw[, i], module_mean)
    if (orient && is.finite(orientation$raw_ME_vs_module_mean_cor[[i]]) &&
        orientation$raw_ME_vs_module_mean_cor[[i]] < 0) {
      orientation$sign_multiplier[[i]] <- -1L
      oriented[, i] <- -raw[, i]
    }
  }

  metadata <- WGCNA_align_sample_info(sample_info, rownames(oriented))
  zscore <- scale(oriented)
  long <- do.call(rbind, lapply(seq_len(ncol(oriented)), function(i) {
    data.frame(
      Sample_ID = rownames(oriented),
      module = module_names[[i]],
      eigengene = oriented[, i],
      eigengene_zscore = zscore[, i],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  if (!is.null(metadata)) {
    long <- cbind(long, metadata[match(long$Sample_ID, rownames(metadata)), , drop = FALSE])
  }

  time_summary <- NULL
  if (!is.null(time_col)) {
    if (is.null(metadata) || !time_col %in% names(metadata)) {
      stop("`time_col` is not present in aligned `sample_info`.", call. = FALSE)
    }
    groups <- unique(metadata[[time_col]])
    time_summary <- do.call(rbind, lapply(seq_len(ncol(oriented)), function(i) {
      do.call(rbind, lapply(groups, function(group) {
        index <- metadata[[time_col]] == group
        data.frame(
          module = module_names[[i]],
          time_group = as.character(group),
          n = sum(index),
          mean_eigengene = mean(oriented[index, i]),
          sd_eigengene = if (sum(index) > 1L) stats::sd(oriented[index, i]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }

  list(
    MEs_raw = raw,
    MEs_oriented = oriented,
    MEs_zscore = zscore,
    orientation = orientation,
    sample_values = long,
    time_summary = time_summary,
    sample_info = metadata
  )
}

#' Calculate module-trait correlations
#'
#' @param modules A `WGCNA_modules` object.
#' @param traits Numeric sample-by-trait data frame with sample row names or a
#'   `Sample_ID` column.
#' @param MEs Optional eigengene matrix, for example oriented eigengenes from
#'   [WGCNA_module_sample()]. Defaults to raw WGCNA eigengenes.
#' @param use Correlation missing-value policy.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#'
#' @return Correlation, P-value, adjusted P-value matrices and a long table.
#' @export
WGCNA_module_trait <- function(
    modules,
    traits,
    MEs = NULL,
    use = "pairwise.complete.obs",
    p_adjust_method = "bonferroni") {

  # Required R packages: WGCNA
  # Base/recommended packages: stats
  WGCNA_require_namespace()
  if (!inherits(modules, "WGCNA_modules")) {
    stop("`modules` must be returned by WGCNA_build_modules().", call. = FALSE)
  }
  if (is.null(MEs)) MEs <- modules$MEs_raw
  traits <- WGCNA_align_sample_info(traits, rownames(MEs))
  if (is.null(traits) || !ncol(traits)) stop("`traits` must not be empty.", call. = FALSE)
  numeric_traits <- vapply(traits, is.numeric, logical(1))
  if (!all(numeric_traits)) {
    stop("Every trait column must be numeric; encode categorical traits explicitly.", call. = FALSE)
  }
  cor_matrix <- stats::cor(MEs, traits, use = use)
  p_matrix <- matrix(NA_real_, nrow(cor_matrix), ncol(cor_matrix), dimnames = dimnames(cor_matrix))
  for (i in seq_len(nrow(cor_matrix))) {
    for (j in seq_len(ncol(cor_matrix))) {
      complete <- stats::complete.cases(MEs[, i], traits[, j])
      n_complete <- sum(complete)
      if (n_complete >= 3L && is.finite(cor_matrix[i, j])) {
        p_matrix[i, j] <- WGCNA::corPvalueStudent(cor_matrix[i, j], n_complete)
      }
    }
  }
  padj_matrix <- matrix(
    stats::p.adjust(as.vector(p_matrix), method = p_adjust_method),
    nrow = nrow(p_matrix), ncol = ncol(p_matrix), dimnames = dimnames(p_matrix)
  )
  long <- do.call(rbind, lapply(seq_len(nrow(cor_matrix)), function(i) {
    data.frame(
      module = rownames(cor_matrix)[[i]],
      trait = colnames(cor_matrix),
      cor = cor_matrix[i, ],
      pvalue = p_matrix[i, ],
      padj = padj_matrix[i, ],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  list(correlation = cor_matrix, pvalue = p_matrix, padj = padj_matrix, table = long)
}

#' Calculate gene-module membership (kME)
#'
#' @param modules A `WGCNA_modules` object.
#' @param MEs Optional eigengene matrix. Defaults to oriented eigengenes.
#' @param p_adjust_method Multiple-testing method.
#'
#' @return Full gene-by-module matrices, a long table, and one assigned-module
#'   record per gene with within-module rank.
#' @export
WGCNA_module_membership <- function(
    modules,
    MEs = NULL,
    p_adjust_method = "bonferroni") {

  # Required R packages: WGCNA
  # Base/recommended packages: stats
  WGCNA_require_namespace()
  if (!inherits(modules, "WGCNA_modules")) {
    stop("`modules` must be returned by WGCNA_build_modules().", call. = FALSE)
  }
  if (is.null(MEs)) MEs <- WGCNA_module_sample(modules)$MEs_oriented
  cor_matrix <- stats::cor(modules$datExpr, MEs, use = "pairwise.complete.obs")
  p_matrix <- WGCNA::corPvalueStudent(cor_matrix, nrow(modules$datExpr))
  padj_matrix <- matrix(
    stats::p.adjust(as.vector(p_matrix), method = p_adjust_method),
    nrow = nrow(p_matrix), ncol = ncol(p_matrix), dimnames = dimnames(p_matrix)
  )
  long <- do.call(rbind, lapply(seq_len(ncol(cor_matrix)), function(i) {
    data.frame(
      GeneID = rownames(cor_matrix),
      ME = colnames(cor_matrix)[[i]],
      module = substring(colnames(cor_matrix)[[i]], 3),
      cor = cor_matrix[, i],
      pvalue = p_matrix[, i],
      padj = padj_matrix[, i],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  assigned <- modules$module_gene
  assigned$ME <- paste0("ME", assigned$module)
  key <- cbind(match(assigned$GeneID, rownames(cor_matrix)), match(assigned$ME, colnames(cor_matrix)))
  assigned$kME <- cor_matrix[key]
  assigned$pvalue <- p_matrix[key]
  assigned$padj <- padj_matrix[key]
  assigned$abs_kME <- abs(assigned$kME)
  assigned$rank_in_module <- ave(
    -assigned$abs_kME, assigned$module,
    FUN = function(x) rank(x, ties.method = "min")
  )
  assigned <- assigned[order(assigned$module, assigned$rank_in_module), , drop = FALSE]
  rownames(assigned) <- NULL
  list(correlation = cor_matrix, pvalue = p_matrix, padj = padj_matrix,
       table = long, assigned = assigned)
}

#' Export intramodular TOM networks
#'
#' Recomputes TOM within each requested module, matching the scalable behavior
#' of the reference image. Grey genes are skipped by default.
#'
#' @param modules A `WGCNA_modules` object.
#' @param module Optional module names; `NULL` exports all eligible modules.
#' @param threshold Minimum TOM weight.
#' @param max_module_genes Maximum module size to export.
#' @param max_edges Optional maximum number of strongest edges per module.
#' @param include_grey Whether to include the grey module.
#' @param output_dir Optional directory for `.nodes.txt` and `.edges.txt` files.
#'
#' @return Named list of node and edge data frames.
#' @export
WGCNA_export_network <- function(
    modules,
    module = NULL,
    threshold = 0.15,
    max_module_genes = 2000,
    max_edges = Inf,
    include_grey = FALSE,
    output_dir = NULL) {

  # Required R packages: WGCNA
  # Base packages: utils
  WGCNA_require_namespace()
  if (!inherits(modules, "WGCNA_modules")) {
    stop("`modules` must be returned by WGCNA_build_modules().", call. = FALSE)
  }
  WGCNA_assert_probability(threshold, "threshold")
  WGCNA_assert_positive(max_module_genes, "max_module_genes")
  WGCNA_assert_flag(include_grey, "include_grey")
  available <- unique(unname(modules$module_colors))
  if (is.null(module)) module <- available
  unknown <- setdiff(module, available)
  if (length(unknown)) stop("Unknown modules: ", paste(unknown, collapse = ", "), call. = FALSE)
  if (!include_grey) module <- setdiff(module, "grey")
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output <- list()
  for (module_name in module) {
    genes <- names(modules$module_colors)[modules$module_colors == module_name]
    if (length(genes) < 2L || length(genes) > max_module_genes) next
    adjacency <- WGCNA::adjacency(
      modules$datExpr[, genes, drop = FALSE],
      power = modules$parameters$power,
      type = modules$parameters$network_type
    )
    TOM <- WGCNA::TOMsimilarity(adjacency, TOMType = modules$parameters$TOM_type)
    indices <- which(upper.tri(TOM) & TOM >= threshold, arr.ind = TRUE)
    edges <- data.frame(
      fromNode = genes[indices[, 1]],
      toNode = genes[indices[, 2]],
      weight = TOM[indices],
      direction = "undirected",
      stringsAsFactors = FALSE
    )
    edges <- edges[order(edges$weight, decreasing = TRUE), , drop = FALSE]
    if (is.finite(max_edges)) edges <- utils::head(edges, as.integer(max_edges))
    nodes <- data.frame(
      nodeName = genes,
      nodeAttr = module_name,
      stringsAsFactors = FALSE
    )
    output[[module_name]] <- list(nodes = nodes, edges = edges)
    if (!is.null(output_dir)) {
      utils::write.table(
        nodes, file.path(output_dir, paste0(module_name, ".nodes.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
      utils::write.table(
        edges, file.path(output_dir, paste0(module_name, ".edges.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }
  output
}

#' Run the reusable WGCNA workflow
#'
#' Orchestrates expression preparation, power selection, module construction,
#' no-trait module-sample summaries, optional module-trait analysis, kME, and
#' optional network export.
#'
#' @param expression Gene-by-sample expression matrix.
#' @param sample_info Optional sample metadata.
#' @param traits Optional numeric sample-by-trait table.
#' @param prepare_args,power_args,module_args Named arguments forwarded to the
#'   corresponding `WGCNA_` functions.
#' @param time_col Optional time/group column in `sample_info`.
#' @param output_dir Optional directory for core CSV/TSV outputs.
#' @param export_networks Whether intramodular networks should be exported.
#' @param network_args Named arguments passed to [WGCNA_export_network()].
#' @param verbose Print a completion summary.
#'
#' @return Invisibly, a `WGCNA_workflow` list of all intermediate and final
#'   computation objects.
#' @export
WGCNA_run <- function(
    expression,
    sample_info = NULL,
    traits = NULL,
    prepare_args = list(),
    power_args = list(),
    module_args = list(),
    time_col = NULL,
    output_dir = NULL,
    export_networks = FALSE,
    network_args = list(),
    verbose = TRUE) {

  # Required R packages: WGCNA
  # This entry point delegates to the other WGCNA_ functions above.
  for (item in c("prepare_args", "power_args", "module_args", "network_args")) {
    if (!is.list(get(item))) stop("`", item, "` must be a list.", call. = FALSE)
  }
  prepared <- do.call(WGCNA_prepare_expression, c(list(expression = expression), prepare_args))
  power <- do.call(WGCNA_select_power, c(list(prepared = prepared), power_args))
  if (!"network_type" %in% names(module_args)) module_args$network_type <- power$network_type
  modules <- do.call(WGCNA_build_modules, c(list(prepared = prepared, power = power), module_args))
  module_sample <- WGCNA_module_sample(modules, sample_info, time_col)
  membership <- WGCNA_module_membership(modules, module_sample$MEs_oriented)
  module_trait <- if (is.null(traits)) NULL else {
    WGCNA_module_trait(modules, traits, MEs = module_sample$MEs_oriented)
  }
  networks <- NULL
  if (export_networks) {
    if (!"output_dir" %in% names(network_args) && !is.null(output_dir)) {
      network_args$output_dir <- file.path(output_dir, "Network")
    }
    networks <- do.call(WGCNA_export_network, c(list(modules = modules), network_args))
  }

  result <- list(
    prepared = prepared,
    power = power,
    modules = modules,
    module_sample = module_sample,
    membership = membership,
    module_trait = module_trait,
    networks = networks,
    call = match.call()
  )
  class(result) <- c("WGCNA_workflow", "list")
  if (!is.null(output_dir)) WGCNA_write_results(result, output_dir)
  if (isTRUE(verbose)) {
    cat("WGCNA_run completed\n")
    cat("  Samples:", nrow(prepared$datExpr), " Genes:", ncol(prepared$datExpr), "\n")
    cat("  Power:", power$power, " Modules including grey:",
        length(unique(modules$module_colors)), "\n")
  }
  invisible(result)
}

#' Write core WGCNA workflow results
#'
#' @param result A `WGCNA_workflow` object.
#' @param output_dir Destination directory.
#'
#' @return Invisibly returns normalized `output_dir`.
#' @export
WGCNA_write_results <- function(result, output_dir) {
  # Required R packages: none
  # Base packages: utils
  if (!inherits(result, "WGCNA_workflow")) {
    stop("`result` must be returned by WGCNA_run().", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(x, filename, row_names = FALSE) {
    utils::write.csv(x, file.path(output_dir, filename), row.names = row_names, quote = FALSE)
  }
  write_csv(result$prepared$qc, "expression_qc.csv")
  write_csv(result$power$fit_indices, "soft_threshold_fit.csv")
  write_csv(result$modules$module_gene, "modulegene.csv")
  write_csv(result$modules$module_summary, "module_summary.csv")
  write_csv(as.data.frame(result$module_sample$MEs_raw), "module_eigengene_raw.csv", TRUE)
  write_csv(as.data.frame(result$module_sample$MEs_oriented), "module_eigengene_oriented.csv", TRUE)
  write_csv(result$module_sample$orientation, "module_eigengene_orientation.csv")
  write_csv(result$module_sample$sample_values, "module_sample_eigengene_long.csv")
  if (!is.null(result$module_sample$time_summary)) {
    write_csv(result$module_sample$time_summary, "module_time_summary.csv")
  }
  write_csv(result$membership$table, "gene_module_correlation.csv")
  write_csv(result$membership$assigned, "assigned_module_kME.csv")
  if (!is.null(result$module_trait)) {
    write_csv(result$module_trait$table, "module_trait_correlation.csv")
  }
  manifest <- data.frame(
    item = c("samples", "genes", "power", "network_type", "modules_including_grey"),
    value = c(
      nrow(result$prepared$datExpr), ncol(result$prepared$datExpr),
      result$power$power, result$power$network_type,
      length(unique(result$modules$module_colors))
    ),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, "RUN_MANIFEST.csv")
  invisible(normalizePath(output_dir, mustWork = TRUE))
}

#' Plot soft-threshold diagnostics
#'
#' @param power_result A `WGCNA_power` object.
#' @param output_file Optional PDF or PNG filename.
#' @param width,height Device dimensions in inches.
#'
#' @return Invisibly returns `power_result`.
#' @export
WGCNA_plot_soft_threshold <- function(
    power_result,
    output_file = NULL,
    width = 16,
    height = 8) {

  # Required R packages: none
  # Base packages: graphics, grDevices, tools
  if (!inherits(power_result, "WGCNA_power")) {
    stop("`power_result` must be returned by WGCNA_select_power().", call. = FALSE)
  }
  if (!is.null(output_file)) WGCNA_open_device(output_file, width, height)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old)
    if (!is.null(output_file)) grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfcol = c(1, 2))
  fit <- power_result$fit_indices
  graphics::plot(
    fit[[1]], power_result$signed_fit, type = "n", ylim = c(0, 1),
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R^2",
    main = "Scale independence"
  )
  graphics::text(fit[[1]], power_result$signed_fit, labels = fit[[1]], col = "red")
  graphics::abline(h = power_result$fit_cutoff, col = "red")
  graphics::abline(v = power_result$power, col = "blue", lty = 2)
  graphics::plot(
    fit[[1]], fit[[5]], type = "n",
    xlab = "Soft Threshold (power)", ylab = "Mean connectivity",
    main = "Mean connectivity"
  )
  graphics::text(fit[[1]], fit[[5]], labels = fit[[1]], col = "red")
  graphics::abline(v = power_result$power, col = "blue", lty = 2)
  invisible(power_result)
}

WGCNA_require_namespace <- function() {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required for this function.", call. = FALSE)
  }
}

WGCNA_assert_flag <- function(value, argument) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", argument, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

WGCNA_assert_probability <- function(value, argument) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value > 1) {
    stop("`", argument, "` must be a number between 0 and 1.", call. = FALSE)
  }
}

WGCNA_assert_positive <- function(value, argument) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value <= 0) {
    stop("`", argument, "` must be one positive number.", call. = FALSE)
  }
}

WGCNA_get_datExpr <- function(x) {
  datExpr <- if (inherits(x, "WGCNA_prepared")) x$datExpr else x
  datExpr <- as.data.frame(datExpr, check.names = FALSE)
  if (!nrow(datExpr) || !ncol(datExpr) || anyNA(datExpr) ||
      any(!is.finite(as.matrix(datExpr)))) {
    stop("A complete finite sample-by-gene matrix is required.", call. = FALSE)
  }
  datExpr
}

WGCNA_image_power <- function(n_samples, network_type) {
  unsigned <- identical(network_type, "unsigned")
  if (n_samples < 20) return(if (unsigned) 9 else 18)
  if (n_samples < 30) return(if (unsigned) 8 else 16)
  if (n_samples < 40) return(if (unsigned) 7 else 14)
  if (unsigned) 6 else 12
}

WGCNA_align_sample_info <- function(sample_info, sample_names) {
  if (is.null(sample_info)) return(NULL)
  info <- as.data.frame(sample_info, check.names = FALSE, stringsAsFactors = FALSE)
  if ("Sample_ID" %in% names(info)) {
    ids <- as.character(info$Sample_ID)
    info$Sample_ID <- NULL
    rownames(info) <- ids
  }
  if (is.null(rownames(info)) || !setequal(rownames(info), sample_names)) {
    stop("Sample metadata IDs must match expression sample names.", call. = FALSE)
  }
  info[sample_names, , drop = FALSE]
}

WGCNA_open_device <- function(filename, width, height) {
  extension <- tolower(tools::file_ext(filename))
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (extension == "pdf") {
    grDevices::cairo_pdf(filename, width = width, height = height)
  } else if (extension == "png") {
    grDevices::png(filename, width = width, height = height, units = "in", res = 300)
  } else {
    stop("`output_file` must end in .pdf or .png.", call. = FALSE)
  }
}
# All WGCNA plot text prefers Times New Roman where a font family is configurable.

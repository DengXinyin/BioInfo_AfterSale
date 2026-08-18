# ============================================================================
# LEfSe analysis and report-style visualisation
#
# Required R packages (install before running these functions):
#   install.packages(c("ggplot2", "dplyr"))
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("microbiomeMarker", "phyloseq", "ggtree"))
#
# To reproduce the output of `192.168.30.202:23099/micro_dy_gro/micro:v2.45`
# exactly, use the same versions of microbiomeMarker, phyloseq, ggtree and
# ggplot2 as installed in that image.  The functions below do NOT call Docker.
# A system-installed font (for example Times New Roman) is additionally needed
# when `font_family = "Times New Roman"` is used for PDF output.
# ============================================================================
#'
#' LEfSe analysis and report-style visualisation
#'
#' Reusable LEfSe functions extracted and refactored from the internal image
#' `192.168.30.202:23099/micro_dy_gro/micro:v2.45` (package
#' `microbiomeMarker`).  This file does not invoke Docker: it runs anywhere R
#' has the packages listed in `DESCRIPTION` installed.
#'
#' Required packages: microbiomeMarker, phyloseq, ggtree, ggplot2 and dplyr.
#' `ggtree` is installed as a dependency of microbiomeMarker but is listed here
#' because cladogram labels are ggtree layers.
#'
#' @param otu Numeric matrix/data.frame with features in rows and samples in
#'   columns. Row names must be feature IDs.
#' @param taxonomy Taxonomy matrix/data.frame with the same feature IDs as
#'   `otu`. Taxonomic ranks should be columns in Kingdom-to-Genus order.
#' @param group Named vector of group labels, or a one-column data.frame whose
#'   row names are sample IDs.
#' @param group_name Name assigned to the grouping variable in phyloseq.
#' @param lda_cutoff Minimum LDA score retained by LEfSe.
#' @param kw_cutoff,wilcoxon_cutoff LEfSe significance cutoffs.
#' @param min_total_abundance Remove features whose total abundance is at or
#'   below this value before LEfSe.
#' @param min_relative_mean Remove features whose mean abundance divided by
#'   the sum of mean abundance is at or below this threshold.
#' @param min_prevalence Keep features present above `min_presence` in more
#'   than `min_prevalence` of samples.
#' @param min_presence Abundance threshold used by `min_prevalence`.
#' @param min_variance Remove features at or below this sample-wise abundance
#'   variance before LEfSe. Set to `NULL` to skip variance filtering.
#'
#' @return A `microbiomeMarker` LEfSe object.
#' @export
LEFSE_run <- function(
    otu,
    taxonomy,
    group,
    group_name = "Group",
    lda_cutoff = 2,
    kw_cutoff = 0.05,
    wilcoxon_cutoff = 0.05,
    min_total_abundance = 10,
    min_relative_mean = 1e-5,
    min_prevalence = 0.2,
    min_presence = 2,
    min_variance = 1e-5) {

  .lefse_require(c("microbiomeMarker", "phyloseq"))
  otu <- .lefse_numeric_matrix(otu, "otu")
  taxonomy <- as.matrix(taxonomy)
  if (is.null(rownames(taxonomy))) {
    stop("`taxonomy` must have feature IDs as row names.", call. = FALSE)
  }
  if (is.null(colnames(taxonomy))) {
    colnames(taxonomy) <- paste0("Rank", seq_len(ncol(taxonomy)))
  }
  common_features <- intersect(rownames(otu), rownames(taxonomy))
  if (!length(common_features)) {
    stop("`otu` and `taxonomy` have no shared feature IDs.", call. = FALSE)
  }
  otu <- otu[common_features, , drop = FALSE]
  taxonomy <- taxonomy[common_features, , drop = FALSE]
  if (!is.null(min_relative_mean)) {
    relative_mean <- rowMeans(otu) / sum(rowMeans(otu))
    keep <- is.finite(relative_mean) & relative_mean > min_relative_mean
    otu <- otu[keep, , drop = FALSE]
    taxonomy <- taxonomy[rownames(otu), , drop = FALSE]
  }
  group <- .lefse_group_data(group, colnames(otu), group_name)

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(group)
  )
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > min_total_abundance, ps)
  if (!is.null(min_variance)) {
    ps <- phyloseq::filter_taxa(ps, function(x) stats::var(x) > min_variance, prune = TRUE)
  }
  ps <- phyloseq::filter_taxa(
    ps,
    function(x) sum(x > min_presence) > min_prevalence * length(x),
    prune = TRUE
  )
  if (!phyloseq::ntaxa(ps)) stop("No features remain after filtering.", call. = FALSE)

  microbiomeMarker::run_lefse(
    ps,
    group = group_name,
    wilcoxon_cutoff = wilcoxon_cutoff,
    kw_cutoff = kw_cutoff,
    multigrp_strat = TRUE,
    lda_cutoff = lda_cutoff
  )
}

#' Extract a portable LEfSe marker table
#'
#' @param lefse_result Result returned by [LEFSE_run()].
#' @return A data frame with feature, enriched group, LDA score and P values.
#' @export
LEFSE_marker_table <- function(lefse_result) {
  .lefse_require("microbiomeMarker")
  as.data.frame(microbiomeMarker::marker_table(lefse_result))
}

#' Plot a LEfSe LDA score bar chart
#'
#' @param marker_table Marker data frame from [LEFSE_marker_table()].
#' @param group_colors Named colours. Names must contain all enriched groups.
#' @param output_file Optional PDF/PNG output path.
#' @param font_family Plot font family.
#' @param lda_cutoff Keep marker rows at or above this LDA score.
#' @param figure_width,figure_height Dimensions in inches. `figure_height =
#'   NULL` scales with number of markers.
#' @return A ggplot object, invisibly saved when `output_file` is supplied.
#' @export
LEFSE_LDA_plot <- function(
    marker_table,
    group_colors,
    output_file = NULL,
    font_family = "Times New Roman",
    lda_cutoff = 2,
    figure_width = 17,
    figure_height = NULL) {

  .lefse_require(c("ggplot2", "dplyr"))
  .lefse_register_font_metrics(font_family)
  marker_table <- .lefse_validate_markers(marker_table)
  d <- marker_table[marker_table$ef_lda >= lda_cutoff, , drop = FALSE]
  if (!nrow(d)) stop("No markers meet `lda_cutoff`.", call. = FALSE)
  group_colors <- .lefse_validate_colors(group_colors, unique(d$enrich_group))
  d$enrich_group <- factor(d$enrich_group, levels = names(group_colors))
  d$Taxonomy <- sub(".*\\|", "", d$feature)
  d <- dplyr::arrange(d, enrich_group, dplyr::desc(ef_lda))
  d$Taxonomy <- factor(
    make.unique(as.character(d$Taxonomy)),
    levels = rev(make.unique(as.character(d$Taxonomy)))
  )
  figure_height <- figure_height %||% max(11, 0.31 * nrow(d) + 4)
  p <- ggplot2::ggplot(d, ggplot2::aes(ef_lda, Taxonomy, fill = enrich_group)) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::scale_fill_manual(values = group_colors, name = NULL, drop = FALSE) +
    ggplot2::geom_vline(
      xintercept = seq(0, ceiling(max(d$ef_lda)), 1),
      linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::labs(title = "LEfSe LDA score", x = "LDA SCORE (log 10)", y = NULL) +
    ggplot2::theme_classic(base_family = font_family, base_size = 20) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 26),
      axis.text.y = ggplot2::element_text(size = 17),
      axis.text.x = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 22),
      legend.position = "top", legend.text = ggplot2::element_text(size = 18),
      plot.margin = ggplot2::margin(14, 38, 14, 14)
    )
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = figure_width, height = figure_height,
                    device = grDevices::cairo_pdf, family = font_family, limitsize = FALSE)
  }
  p
}

#' Plot a LEfSe cladogram
#'
#' @param lefse_result Result returned by [LEFSE_run()].
#' @param group_colors Named colours matching enriched group names.
#' @param output_file Optional PDF output path.
#' @param clade_label_level Taxonomic level at which branch letters begin.
#' @param font_family Font for all text layers.
#' @param label_point_size Size of the coloured taxonomy-label squares.
#' @param figure_width,figure_height Plot dimensions in inches.
#' @return A ggplot/ggtree object, invisibly saved when `output_file` is supplied.
#' @export
LEFSE_cladogram_plot <- function(
    lefse_result,
    group_colors,
    output_file = NULL,
    clade_label_level = 6,
    font_family = "Times New Roman",
    label_point_size = 8,
    figure_width = 28,
    figure_height = 18) {

  .lefse_require(c("microbiomeMarker", "ggplot2", "ggtree"))
  .lefse_register_font_metrics(font_family)
  markers <- LEFSE_marker_table(lefse_result)
  if (!nrow(markers)) stop("`lefse_result` contains no markers.", call. = FALSE)
  group_colors <- .lefse_validate_colors(group_colors, unique(markers$enrich_group))

  p <- microbiomeMarker::plot_cladogram(
    lefse_result,
    color = group_colors,
    clade_label_level = clade_label_level,
    clade_label_font_size = 5.5,
    annotation_shape_size = label_point_size,
    group_legend_param = list(nrow = 1, ncol = length(group_colors), byrow = TRUE),
    marker_legend_param = list(ncol = 2, byrow = TRUE)
  ) +
    ggplot2::theme(
      legend.position = "right", legend.box = "vertical",
      text = ggplot2::element_text(family = font_family),
      legend.text = ggplot2::element_text(size = 19, family = font_family),
      legend.title = ggplot2::element_text(size = 21, family = font_family),
      plot.margin = ggplot2::margin(14, 34, 14, 14)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      nrow = 1, ncol = length(group_colors), byrow = TRUE
    ))

  # microbiomeMarker uses ggtree StatCladeText layers whose default family is
  # sans. Override only these text layers; changing non-text tree layers can
  # break older ggtree releases.
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$stat, "StatCladeText") ||
        inherits(p$layers[[i]]$stat, "StatCladeText2")) {
      p$layers[[i]]$aes_params$family <- font_family
    }
  }
  # The generic one-row guide is needed by older ggplot2 builds. Restore the
  # exact colours because that guide otherwise falls back to the default palette.
  p$guides$guides$fill$params$override.aes$fill <- unname(group_colors)

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = figure_width, height = figure_height,
                    device = grDevices::cairo_pdf, limitsize = FALSE)
  }
  p
}

.lefse_require <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    stop("Install required package(s): ", paste(absent, collapse = ", "), call. = FALSE)
  }
}

.lefse_register_font_metrics <- function(font_family) {
  # `cairo_pdf()` resolves the actual font through fontconfig, while grid needs
  # a known PDF metric family before a device is opened.  This alias lets a
  # system-installed Times New Roman be used on standalone R installations.
  aliases <- c("Times New Roman" = "Times", "Arial" = "Helvetica")
  if (!font_family %in% names(aliases) || font_family %in% names(grDevices::pdfFonts())) {
    return(invisible(NULL))
  }
  metric <- grDevices::pdfFonts(unname(aliases[[font_family]]))[[1L]]
  do.call(grDevices::pdfFonts, stats::setNames(list(metric), font_family))
  invisible(NULL)
}

.lefse_numeric_matrix <- function(x, argument) {
  x <- as.matrix(x)
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("`", argument, "` needs row and column names.", call. = FALSE)
  }
  suppressWarnings(storage.mode(x) <- "double")
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("`", argument, "` must contain finite non-negative abundances.", call. = FALSE)
  }
  x
}

.lefse_group_data <- function(group, sample_ids, group_name) {
  if (is.data.frame(group)) {
    if (ncol(group) != 1L) stop("`group` data.frame must have one column.", call. = FALSE)
    group_ids <- rownames(group)
    group <- group[[1L]]
    names(group) <- group_ids
  }
  if (is.null(names(group))) {
    if (length(group) != length(sample_ids)) {
      stop("Unnamed `group` must have one value per sample.", call. = FALSE)
    }
    names(group) <- sample_ids
  }
  if (length(missing <- setdiff(sample_ids, names(group)))) {
    stop("`group` is missing sample(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(stringsAsFactors = FALSE, row.names = sample_ids)
  out[[group_name]] <- factor(as.character(group[sample_ids]))
  if (anyNA(out[[group_name]]) || any(!nzchar(as.character(out[[group_name]])))) {
    stop("`group` cannot contain missing or empty labels.", call. = FALSE)
  }
  out
}

.lefse_validate_markers <- function(marker_table) {
  marker_table <- as.data.frame(marker_table, stringsAsFactors = FALSE)
  needed <- c("feature", "enrich_group", "ef_lda")
  if (length(missing <- setdiff(needed, names(marker_table)))) {
    stop("`marker_table` is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  marker_table$ef_lda <- as.numeric(marker_table$ef_lda)
  marker_table
}

.lefse_validate_colors <- function(group_colors, groups) {
  if (!is.character(group_colors) || is.null(names(group_colors))) {
    stop("`group_colors` must be a named character vector.", call. = FALSE)
  }
  missing <- setdiff(as.character(groups), names(group_colors))
  if (length(missing)) {
    stop("`group_colors` is missing group(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  group_colors[unique(as.character(groups))]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

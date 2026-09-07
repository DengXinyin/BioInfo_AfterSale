# Differential abundance visualizations. Existing result tables can be passed directly.

metage_diff_boxplot <- function(abundance, group, feature, p_value = NULL,
                                 method = c("none", "wilcox", "anova"), relative = TRUE,
                                 colors = NULL, font_family = metage_font_family()) {
  method <- match.arg(method)
  x <- .metage_matrix(abundance)
  if (!feature %in% rownames(x)) stop("Unknown feature: ", feature)
  groups <- .metage_group(group, colnames(x)); value <- x[feature, ]
  if (relative) value <- value / pmax(colSums(x), 1)
  dat <- data.frame(Sample = colnames(x), Group = groups, Value = value)
  if (is.null(p_value) && method == "wilcox" && nlevels(groups) == 2)
    p_value <- stats::wilcox.test(Value ~ Group, dat)$p.value
  if (is.null(p_value) && method == "anova" && nlevels(groups) > 1)
    p_value <- summary(stats::aov(Value ~ Group, dat))[[1]][["Pr(>F)"]][1]
  if (is.null(colors)) colors <- .metage_palette(nlevels(groups))
  ggplot2::ggplot(dat, ggplot2::aes(Group, Value, fill = Group)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, size = 2) + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = if (relative) "Relative abundance" else "Abundance",
                  title = feature,
                  subtitle = if (!is.null(p_value)) paste0("P = ", format.pval(p_value, 3)) else NULL) +
    .metage_theme(font_family) + ggplot2::theme(legend.position = "none")
}

metage_stamp_plot <- function(result, feature = "Feature", difference = "Difference",
                               lower = "Lower", upper = "Upper", p = "P",
                               p_adjusted = FALSE, top_n = 20,
                               font_family = metage_font_family()) {
  required <- c(feature, difference, lower, upper, p)
  if (!all(required %in% names(result))) stop("STAMP result columns are missing.")
  dat <- result
  for (column in c(difference, lower, upper, p))
    dat[[column]] <- suppressWarnings(as.numeric(dat[[column]]))
  dat <- dat[stats::complete.cases(dat[required]), , drop = FALSE]
  dat <- dat[order(dat[[p]]), , drop = FALSE]
  dat <- utils::head(dat, top_n)
  dat$.Feature <- factor(dat[[feature]], levels = rev(dat[[feature]]))
  dat$.Difference <- as.numeric(dat[[difference]])
  dat$.Lower <- as.numeric(dat[[lower]]); dat$.Upper <- as.numeric(dat[[upper]])
  ggplot2::ggplot(dat, ggplot2::aes(.Feature, .Difference)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .Lower, ymax = .Upper), width = 0.25) +
    ggplot2::geom_point(size = 2.6, colour = "#377EB8") + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Difference between proportions",
                  title = "STAMP effect size",
                  subtitle = paste(if (p_adjusted) "Adjusted" else "Raw", "P-value column:", p)) +
    .metage_theme(font_family)
}

metage_random_forest_importance <- function(importance, feature = "Feature",
                                             value = "MeanDecreaseGini", top_n = 20,
                                             font_family = metage_font_family()) {
  if (!all(c(feature, value) %in% names(importance))) stop("Importance columns are missing.")
  importance[[value]] <- suppressWarnings(as.numeric(importance[[value]]))
  importance <- importance[is.finite(importance[[value]]), , drop = FALSE]
  dat <- importance[order(importance[[value]], decreasing = TRUE), , drop = FALSE]
  dat <- utils::head(dat, top_n); dat$.Feature <- factor(dat[[feature]], levels = rev(dat[[feature]]))
  ggplot2::ggplot(dat, ggplot2::aes(.Feature, .data[[value]])) +
    ggplot2::geom_col(fill = "#377EB8", width = 0.72) + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = value, title = "Random forest feature importance") + .metage_theme(font_family)
}

metage_lda_score_plot <- function(marker, feature = "Feature", group = "Group", score = "LDA",
                                   cutoff = 2, colors = NULL,
                                   font_family = metage_font_family()) {
  if (!all(c(feature, group, score) %in% names(marker))) stop("LEfSe marker columns are missing.")
  marker[[score]] <- suppressWarnings(as.numeric(marker[[score]]))
  dat <- marker[is.finite(marker[[score]]) & abs(marker[[score]]) >= cutoff, , drop = FALSE]
  if (!nrow(dat)) stop("No marker passes `cutoff`.")
  dat$.Feature <- factor(dat[[feature]], levels = dat[[feature]][order(abs(dat[[score]]))])
  dat$.Group <- factor(dat[[group]], levels = unique(dat[[group]])); dat$.Score <- dat[[score]]
  if (is.null(colors)) colors <- .metage_palette(nlevels(dat$.Group))
  ggplot2::ggplot(dat, ggplot2::aes(.Feature, .Score, fill = .Group)) +
    ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "LDA score", title = paste0("LEfSe markers (|LDA| >= ", cutoff, ")")) +
    .metage_theme(font_family)
}

metage_significant_heatmap <- function(abundance, result, feature = "Feature", p = "Padj",
                                        cutoff = 0.05, top_n = 30, ...) {
  if (!all(c(feature, p) %in% names(result))) stop("Result feature/P columns are missing.")
  result[[p]] <- suppressWarnings(as.numeric(result[[p]]))
  keep <- as.character(result[[feature]][is.finite(result[[p]]) & result[[p]] <= cutoff])
  keep <- intersect(keep, rownames(abundance)); keep <- utils::head(keep, top_n)
  if (!length(keep)) stop("No significant feature is present in `abundance`.")
  metage_abundance_heatmap(abundance[keep, , drop = FALSE], top_n = length(keep), ...)
}

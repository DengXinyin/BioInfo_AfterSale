# Differential-feature visualizations consolidated from anova.r,
# plot_wilcoxon.r, tax/func2/fungi_stamp.R, tax/func2/fungi_randomForest.R,
# metagenomeSeq.r, plot_metastats.r, and tax_LDAscore.R in micro_v2:v5.42.

micro_diff_boxplot <- function(abundance, group, feature, method = c("wilcox", "anova", "none"),
                               colors = NULL, relative = TRUE, font_family = "sans") {
  method <- match.arg(method)
  x <- .micro_matrix(abundance)
  if (!feature %in% rownames(x)) stop("Unknown feature: ", feature)
  groups <- .micro_group(group, colnames(x))
  value <- x[feature, ]
  if (relative) value <- value / pmax(colSums(x), 1)
  dat <- data.frame(Sample = colnames(x), Group = groups, Value = value)
  pvalue <- NA_real_
  if (method == "wilcox" && nlevels(groups) == 2) pvalue <- stats::wilcox.test(Value ~ Group, dat)$p.value
  if (method == "anova" && nlevels(groups) > 1) pvalue <- summary(stats::aov(Value ~ Group, dat))[[1]][["Pr(>F)"]][1]
  if (is.null(colors)) colors <- .micro_palette(nlevels(groups))
  ggplot2::ggplot(dat, ggplot2::aes(x = Group, y = Value, fill = Group)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, size = 2) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = if (relative) "Relative abundance" else "Abundance",
                  title = feature,
                  subtitle = if (is.finite(pvalue)) paste0(method, " P = ", format.pval(pvalue, 3)) else NULL) +
    .micro_theme(font_family) + ggplot2::theme(legend.position = "none")
}

micro_stamp_plot <- function(abundance, group, adjust = "BH", top_n = 20,
                             colors = c("#E69F00", "#56B4E9"),
                             font_family = "sans") {
  x <- .micro_matrix(abundance)
  groups <- .micro_group(group, colnames(x))
  if (nlevels(groups) != 2) stop("STAMP-style comparison requires exactly two groups.")
  x <- sweep(x, 2, pmax(colSums(x), 1), "/")
  lev <- levels(groups)
  stats <- lapply(seq_len(nrow(x)), function(i) {
    z <- stats::t.test(x[i, groups == lev[1]], x[i, groups == lev[2]])
    data.frame(Feature = rownames(x)[i], Difference = unname(diff(rev(z$estimate))),
               Lower = z$conf.int[1], Upper = z$conf.int[2], P = z$p.value)
  })
  stats <- do.call(rbind, stats)
  stats$Padj <- stats::p.adjust(stats$P, method = adjust)
  stats <- stats[order(stats$Padj), , drop = FALSE]
  stats <- utils::head(stats, top_n)
  means <- aggregate(cbind(Value = as.vector(x[stats$Feature, , drop = FALSE])),
                     list(Feature = rep(stats$Feature, ncol(x)),
                          Group = rep(groups, each = nrow(stats))), mean)
  means$Feature <- factor(means$Feature, levels = rev(stats$Feature))
  stats$Feature <- factor(stats$Feature, levels = levels(means$Feature))
  p1 <- ggplot2::ggplot(means, ggplot2::aes(x = Feature, y = Value, fill = Group)) +
    ggplot2::geom_col(position = "dodge", width = 0.7, colour = "black", linewidth = 0.2) +
    ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "Mean proportion") + .micro_theme(font_family)
  p2 <- ggplot2::ggplot(stats, ggplot2::aes(x = Feature, y = Difference)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.3) +
    ggplot2::geom_point(size = 2.5, colour = "#2C7FB8") + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = paste(lev, collapse = " minus "),
                  title = "95% confidence intervals") + .micro_theme(font_family) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
  grob <- gridExtra::arrangeGrob(p1, p2, ncol = 2, widths = c(1.2, 1))
  structure(list(plot = grob, statistics = stats, means = means), class = "micro_stamp")
}

micro_random_forest_plot <- function(abundance, group, top_n = 10, seed = 123,
                                     font_family = "sans") {
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("Package `randomForest` is required.")
  x <- .micro_matrix(abundance)
  groups <- .micro_group(group, colnames(x))
  keep <- apply(x, 1, stats::var) > 0
  set.seed(seed)
  fit <- randomForest::randomForest(x = t(x[keep, , drop = FALSE]), y = groups, importance = TRUE)
  imp <- as.data.frame(randomForest::importance(fit))
  metric <- if ("MeanDecreaseGini" %in% names(imp)) "MeanDecreaseGini" else names(imp)[ncol(imp)]
  imp$Feature <- rownames(imp)
  imp <- imp[order(imp[[metric]], decreasing = TRUE), , drop = FALSE]
  imp <- utils::head(imp, top_n)
  imp$Feature <- factor(imp$Feature, levels = rev(imp$Feature))
  imp$Importance <- imp[[metric]]
  p <- ggplot2::ggplot(imp, ggplot2::aes(x = Feature, y = Importance)) +
    ggplot2::geom_col(fill = "#377EB8", width = 0.72) + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, title = "Random forest feature importance") + .micro_theme(font_family)
  structure(list(plot = p, importance = imp, model = fit), class = "micro_random_forest")
}

micro_lda_score_plot <- function(marker, feature = "Feature", group = "Group", score = "LDA",
                                 cutoff = 2, colors = NULL, font_family = "sans") {
  if (!all(c(feature, group, score) %in% names(marker))) stop("Marker columns are missing.")
  dat <- marker[is.finite(marker[[score]]) & abs(marker[[score]]) >= cutoff, , drop = FALSE]
  dat$.Group <- factor(dat[[group]], levels = unique(dat[[group]]))
  dat$.Feature <- factor(dat[[feature]], levels = dat[[feature]][order(abs(dat[[score]]))])
  dat$.Score <- dat[[score]]
  if (is.null(colors)) colors <- .micro_palette(nlevels(dat$.Group))
  ggplot2::ggplot(dat, ggplot2::aes(x = .Feature, y = .Score, fill = .Group)) +
    ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "LDA score", title = paste0("LEfSe markers (LDA >= ", cutoff, ")")) +
    .micro_theme(font_family)
}

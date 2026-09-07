# Alpha diversity and PCA/PCoA/NMDS reconstructed from v2.87.

metage_alpha_indices <- function(abundance) {
  .metage_require("vegan")
  x <- .metage_matrix(abundance)
  observed <- colSums(x > 0)
  # vegan treats rows as sites for matrix input; transpose our features x samples contract.
  est <- vegan::estimateR(t(x))
  data.frame(
    Sample = colnames(x), Observed = observed,
    Chao1 = unname(est["S.chao1", ]), ACE = unname(est["S.ACE", ]),
    Shannon = vegan::diversity(t(x), "shannon"),
    Gini_Simpson = vegan::diversity(t(x), "simpson"),
    check.names = FALSE
  )
}

metage_alpha_boxplot <- function(data, index, group, sample = NULL,
                                  test = c("auto", "wilcox", "anova", "none"),
                                  colors = NULL, font_family = metage_font_family()) {
  test <- match.arg(test)
  if (!all(c(index, group) %in% names(data))) stop("Index/group columns are missing.")
  dat <- data[stats::complete.cases(data[c(index, group)]), , drop = FALSE]
  dat$.Index <- dat[[index]]
  dat$.Group <- factor(dat[[group]], levels = unique(dat[[group]]))
  if (test == "auto") test <- if (nlevels(dat$.Group) == 2) "wilcox" else "anova"
  pvalue <- NA_real_
  if (test == "wilcox" && nlevels(dat$.Group) == 2) {
    pvalue <- stats::wilcox.test(.Index ~ .Group, dat)$p.value
  } else if (test == "anova" && nlevels(dat$.Group) > 1) {
    pvalue <- summary(stats::aov(.Index ~ .Group, dat))[[1]][["Pr(>F)"]][1]
  }
  if (is.null(colors)) colors <- .metage_palette(nlevels(dat$.Group))
  p <- ggplot2::ggplot(dat, ggplot2::aes(.Group, .Index, colour = .Group)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, size = 2.5) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::labs(x = NULL, y = index, title = paste(index, "diversity"),
                  subtitle = if (is.finite(pvalue)) paste0(toupper(test), " P = ", format.pval(pvalue, 3)) else NULL) +
    .metage_theme(font_family) + ggplot2::theme(legend.position = "none")
  if (!is.null(sample) && sample %in% names(dat)) {
    if (requireNamespace("ggrepel", quietly = TRUE))
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[sample]]), size = 3)
  }
  structure(list(plot = p, p_value = pvalue, method = test, data = dat), class = "metage_alpha")
}

metage_ordination <- function(abundance = NULL, group, method = c("PCoA", "NMDS", "PCA"),
                               distance = "bray", coordinates = NULL, labels = FALSE,
                               ellipse = FALSE, colors = NULL, seed = 123,
                               font_family = metage_font_family()) {
  method <- match.arg(method)
  if (!is.null(coordinates)) {
    scores <- as.data.frame(coordinates)
    if (ncol(scores) < 2) stop("`coordinates` needs at least two columns.")
    scores <- scores[, 1:2, drop = FALSE]
    explained <- attr(coordinates, "explained") %metage_or% c(NA_real_, NA_real_)
    fit <- NULL
  } else {
    x <- .metage_matrix(abundance)
    set.seed(seed)
    if (method == "PCA") {
      keep <- apply(x, 1, stats::sd) > 0
      if (sum(keep) < 2) stop("PCA requires at least two variable features.")
      fit <- stats::prcomp(t(x[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)
      scores <- as.data.frame(fit$x[, 1:2, drop = FALSE])
      explained <- 100 * fit$sdev[1:2]^2 / sum(fit$sdev^2)
    } else {
      .metage_require("vegan")
      d <- vegan::vegdist(t(x), method = distance)
      if (method == "PCoA") {
        fit <- stats::cmdscale(d, k = 2, eig = TRUE, add = TRUE)
        scores <- as.data.frame(fit$points)
        explained <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
      } else {
        fit <- vegan::metaMDS(d, k = 2, trymax = 100, trace = FALSE)
        scores <- as.data.frame(fit$points); explained <- c(NA_real_, NA_real_)
      }
    }
  }
  names(scores)[1:2] <- c("Axis1", "Axis2")
  scores$Sample <- rownames(scores) %metage_or% paste0("Sample_", seq_len(nrow(scores)))
  groups <- .metage_group(group, scores$Sample)
  scores$Group <- groups[scores$Sample]
  if (is.null(colors)) colors <- .metage_palette(nlevels(groups))
  axis_label <- function(i) if (is.finite(explained[i])) sprintf("%s%d (%.1f%%)", method, i, explained[i]) else paste0(method, i)
  subtitle <- if (method == "NMDS" && !is.null(fit)) sprintf("Stress = %.3f", fit$stress) else NULL
  p <- ggplot2::ggplot(scores, ggplot2::aes(Axis1, Axis2, colour = Group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey75") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey75") +
    ggplot2::geom_point(size = 3) + ggplot2::scale_colour_manual(values = colors) +
    ggplot2::labs(x = axis_label(1), y = axis_label(2), title = method, subtitle = subtitle) +
    .metage_theme(font_family)
  if (labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Sample))
    else p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = -0.5)
  }
  if (ellipse) {
    if (requireNamespace("ggforce", quietly = TRUE))
      p <- p + ggforce::geom_mark_ellipse(ggplot2::aes(fill = Group), alpha = 0.08, show.legend = FALSE)
    else if (all(table(scores$Group) >= 3)) p <- p + ggplot2::stat_ellipse()
  }
  structure(list(plot = p, scores = scores, model = fit, explained = explained), class = "metage_ordination")
}

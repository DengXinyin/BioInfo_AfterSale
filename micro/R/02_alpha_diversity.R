# Alpha diversity and sampling-depth plots consolidated from a_diversity.r,
# rarefaction.r, Rank_Abundance.r, and Specaccum.r in micro_v2:v5.42.

micro_alpha_boxplot <- function(data, index, group, sample = NULL,
                                test = c("auto", "wilcox", "anova", "none"),
                                colors = NULL, title = NULL, font_family = "sans") {
  test <- match.arg(test)
  if (!all(c(index, group) %in% names(data))) stop("`index` and `group` columns are required.")
  dat <- data[stats::complete.cases(data[, c(index, group)]), , drop = FALSE]
  dat$.Group <- factor(dat[[group]], levels = unique(dat[[group]]))
  dat$.Index <- dat[[index]]
  if (!is.null(sample) && sample %in% names(dat)) dat$.Label <- dat[[sample]]
  if (!is.numeric(dat[[index]])) stop("The diversity index must be numeric.")
  if (is.null(colors)) colors <- .micro_palette(nlevels(dat$.Group))
  if (test == "auto") test <- if (nlevels(dat$.Group) == 2) "wilcox" else "anova"
  subtitle <- NULL
  if (test == "wilcox" && nlevels(dat$.Group) == 2) {
    p <- stats::wilcox.test(.Index ~ .Group, data = dat)$p.value
    subtitle <- paste0("Wilcoxon P = ", format.pval(p, digits = 3))
  } else if (test == "anova" && nlevels(dat$.Group) > 1) {
    fit <- stats::aov(.Index ~ .Group, data = dat)
    p <- summary(fit)[[1]][["Pr(>F)"]][1]
    subtitle <- paste0("ANOVA P = ", format.pval(p, digits = 3))
  }
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .Group, y = .Index, fill = .Group)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.12, size = 2, alpha = 0.75) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = index, title = title %||% paste(index, "diversity"),
                  subtitle = subtitle) + .micro_theme(font_family) +
    ggplot2::theme(legend.position = "none")
  if (!is.null(sample) && sample %in% names(dat) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .Label), size = 3,
                                      max.overlaps = 20, show.legend = FALSE)
  }
  p
}

micro_rarefaction_curve <- function(data, depth, richness, sample, group = NULL,
                                    colors = NULL, title = "Rarefaction curves",
                                    font_family = "sans") {
  required <- c(depth, richness, sample)
  if (!all(required %in% names(data))) stop("Missing rarefaction columns.")
  dat <- data[stats::complete.cases(data[, required]), , drop = FALSE]
  color_col <- if (!is.null(group) && group %in% names(dat)) group else sample
  dat$.Depth <- dat[[depth]]
  dat$.Richness <- dat[[richness]]
  dat$.Sample <- dat[[sample]]
  dat$.Color <- factor(dat[[color_col]], levels = unique(dat[[color_col]]))
  if (is.null(colors)) colors <- .micro_palette(nlevels(dat$.Color))
  ggplot2::ggplot(dat, ggplot2::aes(x = .Depth, y = .Richness, group = .Sample,
                                    colour = .Color)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::labs(x = "Sequencing depth", y = "Observed features", title = title) +
    .micro_theme(font_family)
}

micro_rank_abundance <- function(abundance, group = NULL, colors = NULL,
                                 title = "Rank-abundance curve", font_family = "sans") {
  x <- .micro_matrix(abundance)
  x <- sweep(x, 2, pmax(colSums(x), 1), "/")
  dat <- do.call(rbind, lapply(seq_len(ncol(x)), function(i) {
    y <- sort(x[, i], decreasing = TRUE)
    data.frame(Sample = colnames(x)[i], Rank = seq_along(y), Abundance = y)
  }))
  groups <- .micro_group(group, colnames(x))
  dat$Group <- groups[match(dat$Sample, names(groups))]
  if (is.null(colors)) colors <- .micro_palette(nlevels(dat$Group))
  ggplot2::ggplot(dat, ggplot2::aes(x = Rank, y = Abundance, group = Sample, colour = Group)) +
    ggplot2::geom_line(alpha = 0.65, linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(y = "Relative abundance (log10)", title = title) + .micro_theme(font_family)
}

micro_species_accumulation <- function(abundance, permutations = 100,
                                       title = "Species accumulation curve",
                                       font_family = "sans") {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Package `vegan` is required.")
  x <- .micro_matrix(abundance)
  fit <- vegan::specaccum(t(x), method = "random", permutations = permutations)
  dat <- data.frame(Samples = fit$sites, Richness = fit$richness, SD = fit$sd)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Samples, y = Richness)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Richness - SD, ymax = Richness + SD),
                         fill = "#9ECAE1", alpha = 0.45) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 1) +
    ggplot2::geom_point(colour = "#2166AC", size = 2) +
    ggplot2::labs(y = "Accumulated features", title = title) + .micro_theme(font_family)
  structure(list(plot = p, model = fit, data = dat), class = "micro_accumulation")
}

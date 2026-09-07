# Beta-diversity tests consolidated from Anosim.r, Adonis.r, MRPP.r and the
# tax/func2/fungi variants in micro_v2:v5.42.

micro_beta_test <- function(distance, group, method = c("anosim", "adonis", "mrpp"),
                            permutations = 999, seed = 123) {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Package `vegan` is required.")
  method <- match.arg(method)
  d <- if (inherits(distance, "dist")) distance else stats::as.dist(.micro_matrix(distance, "distance"))
  labels <- attr(d, "Labels")
  groups <- .micro_group(group, labels)
  set.seed(seed)
  if (method == "anosim") vegan::anosim(d, groups, permutations = permutations)
  else if (method == "adonis") vegan::adonis2(d ~ groups, permutations = permutations)
  else vegan::mrpp(d, groups, permutations = permutations)
}

micro_anosim_plot <- function(distance, group, permutations = 999, seed = 123,
                              colors = NULL, font_family = "sans") {
  d <- if (inherits(distance, "dist")) distance else stats::as.dist(.micro_matrix(distance, "distance"))
  labels <- attr(d, "Labels")
  groups <- .micro_group(group, labels)
  fit <- micro_beta_test(d, groups, "anosim", permutations, seed)
  m <- as.matrix(d)
  pair <- which(upper.tri(m), arr.ind = TRUE)
  type <- ifelse(groups[pair[, 1]] == groups[pair[, 2]], "Within groups", "Between groups")
  dat <- data.frame(Type = factor(type, c("Within groups", "Between groups")),
                    Rank = rank(m[pair]))
  if (is.null(colors)) colors <- c("#56B4E9", "#E69F00")
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Type, y = Rank, fill = Type)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.5) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "Ranked dissimilarity", title = "ANOSIM",
                  subtitle = paste0("R = ", round(fit$statistic, 3),
                                    ", P = ", format.pval(fit$signif, 3))) +
    .micro_theme(font_family) + ggplot2::theme(legend.position = "none")
  structure(list(plot = p, model = fit, ranks = dat), class = "micro_anosim")
}

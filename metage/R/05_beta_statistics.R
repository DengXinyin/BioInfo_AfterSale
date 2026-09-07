# Beta-diversity hypothesis tests and ANOSIM visualization.

metage_beta_test <- function(distance, group, method = c("anosim", "adonis", "mrpp"),
                              permutations = 999, seed = 123) {
  .metage_require("vegan"); method <- match.arg(method)
  d <- if (inherits(distance, "dist")) distance else stats::as.dist(.metage_matrix(distance, "distance"))
  groups <- .metage_group(group, attr(d, "Labels"))
  set.seed(seed)
  if (method == "anosim") vegan::anosim(d, groups, permutations = permutations)
  else if (method == "adonis") vegan::adonis2(d ~ groups, permutations = permutations)
  else vegan::mrpp(d, groups, permutations = permutations)
}

metage_anosim_plot <- function(distance, group, model = NULL, permutations = 999,
                                colors = c("#56B4E9", "#E69F00"), seed = 123,
                                font_family = metage_font_family()) {
  d <- if (inherits(distance, "dist")) distance else stats::as.dist(.metage_matrix(distance, "distance"))
  groups <- .metage_group(group, attr(d, "Labels"))
  fit <- model %metage_or% metage_beta_test(d, groups, "anosim", permutations, seed)
  z <- as.matrix(d); pair <- which(upper.tri(z), arr.ind = TRUE)
  dat <- data.frame(Type = ifelse(groups[pair[, 1]] == groups[pair[, 2]],
                                  "Within groups", "Between groups"),
                    Rank = rank(z[pair]))
  dat$Type <- factor(dat$Type, c("Within groups", "Between groups"))
  p <- ggplot2::ggplot(dat, ggplot2::aes(Type, Rank, fill = Type)) +
    ggplot2::geom_boxplot(width = 0.65) + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "Ranked dissimilarity", title = "ANOSIM",
                  subtitle = paste0("R = ", round(fit$statistic, 3),
                                    ", P = ", format.pval(fit$signif, 3))) +
    .metage_theme(font_family) + ggplot2::theme(legend.position = "none")
  structure(list(plot = p, model = fit, ranks = dat), class = "metage_anosim")
}

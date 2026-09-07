#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "micro/tests/test_functions.R"
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
grDevices::pdf(NULL)
source(file.path(root, "micro", "R", "utils.R"))
for (f in sort(list.files(file.path(root, "micro", "R"), pattern = "^[0-9].*[.]R$", full.names = TRUE))) source(f)

set.seed(42)
x <- matrix(rpois(20 * 8, 20), nrow = 20,
            dimnames = list(paste0("F", 1:20), paste0("S", 1:8)))
x[sample(length(x), 30)] <- 0
g <- setNames(rep(c("A", "B"), each = 4), colnames(x))

stopifnot(inherits(micro_abundance_bar(x, g), "ggplot"))
stopifnot(inherits(micro_abundance_heatmap(x), "ggplot"))
stopifnot(inherits(micro_clustered_abundance_bar(x, group = g), "ggplot"))
stopifnot(inherits(micro_rank_abundance(x, g), "ggplot"))
alpha <- data.frame(Sample = colnames(x), Group = unname(g), Shannon = vegan::diversity(t(x)))
stopifnot(inherits(micro_alpha_boxplot(alpha, "Shannon", "Group"), "ggplot"))
curve <- do.call(rbind, lapply(colnames(x), function(s) {
  data.frame(Sample = s, Group = unname(g[s]), Depth = seq(100, 1000, 100),
             Richness = log1p(seq(100, 1000, 100)) * runif(1, 2, 4))
}))
stopifnot(inherits(micro_rarefaction_curve(curve, "Depth", "Richness", "Sample", "Group"), "ggplot"))
stopifnot(inherits(micro_species_accumulation(x, permutations = 10)$plot, "ggplot"))
stopifnot(inherits(micro_ordination(x, g, "PCA")$plot, "ggplot"))
stopifnot(inherits(micro_ordination(x, g, "PCoA")$plot, "ggplot"))
stopifnot(inherits(micro_ordination(x, g, "NMDS")$plot, "ggplot"))
stopifnot(inherits(micro_distance_heatmap(as.matrix(vegan::vegdist(t(x)))), "ggplot"))
stopifnot(inherits(micro_diff_boxplot(x, g, rownames(x)[1]), "ggplot"))
stopifnot(inherits(micro_stamp_plot(x, g)$plot, c("grob", "gtable", "gTree")))
stopifnot(inherits(micro_random_forest_plot(x, g, top_n = 5)$plot, "ggplot"))
markers <- data.frame(Feature = rownames(x)[1:6], Group = rep(c("A", "B"), 3), LDA = seq(2, 4.5, .5))
stopifnot(inherits(micro_lda_score_plot(markers), "ggplot"))
d <- vegan::vegdist(t(x))
stopifnot(inherits(micro_anosim_plot(d, g, permutations = 19)$plot, "ggplot"))
stopifnot(inherits(micro_beta_test(d, g, "adonis", permutations = 19), "anova"))
stopifnot(inherits(micro_overlap_plot(list(A = letters[1:8], B = letters[5:12])), "ggplot"))
edges <- data.frame(Source = paste0("N", 1:8), Target = paste0("N", c(2:8, 1)),
                    Weight = c(.8, -.6, .5, .7, -.4, .9, .3, .6))
stopifnot(inherits(micro_network_plot(edges), "micro_network"))
stopifnot(inherits(micro_zipi_plot(edges)$plot, "ggplot"))
stopifnot(inherits(micro_length_distribution(rnorm(100, 420, 25)), "ggplot"))
stopifnot(inherits(micro_ncm_plot(x), "ggplot"))
grDevices::dev.off()
message("All micro visualization smoke tests passed.")

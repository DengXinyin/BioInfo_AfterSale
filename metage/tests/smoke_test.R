args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
invisible(lapply(list.files(file.path(root, "metage", "R"), pattern = "[.]R$", full.names = TRUE), source))

set.seed(287)
x <- matrix(rpois(18 * 8, 20), 18, 8,
            dimnames = list(paste0("F", 1:18), paste0("S", 1:8)))
g <- setNames(rep(c("A", "B"), each = 4), colnames(x))

stopifnot(inherits(metage_abundance_bar(x), "ggplot"))
stopifnot(inherits(metage_abundance_heatmap(x), "ggplot"))
stopifnot(inherits(metage_sample_correlation(x), "ggplot"))
stopifnot(inherits(metage_read_composition(x[1:3, ]), "ggplot"))
stopifnot(inherits(metage_error_profile(data.frame(Position = 1:10, ErrorRate = runif(10))), "ggplot"))
stopifnot(inherits(metage_length_distribution(sample(500:5000, 100)), "ggplot"))
stopifnot(nrow(metage_alpha_indices(x)) == ncol(x))
stopifnot(inherits(metage_alpha_boxplot(transform(metage_alpha_indices(x), Group = g),
                                       "Shannon", "Group")$plot, "ggplot"))
stopifnot(inherits(metage_ordination(x, g, "PCA")$plot, "ggplot"))
stopifnot(inherits(metage_diff_boxplot(x, g, "F1", method = "wilcox"), "ggplot"))
stopifnot(inherits(metage_beta_test(vegan::vegdist(t(x)), g, "anosim", 9), "anosim"))
stopifnot(inherits(metage_anosim_plot(vegan::vegdist(t(x)), g, permutations = 9)$plot, "ggplot"))

marker <- data.frame(Feature = paste0("F", 1:5), Group = rep(c("A", "B"), length.out = 5),
                     LDA = c(4, 3, -3, 2.5, 1))
stopifnot(inherits(metage_lda_score_plot(marker), "ggplot"))
stamp <- data.frame(Feature = paste0("F", 1:5), Difference = seq(-0.2, 0.2, length.out = 5),
                    Lower = seq(-0.3, 0.1, length.out = 5), Upper = seq(-0.1, 0.3, length.out = 5),
                    P = c(0.001, 0.01, 0.03, 0.2, 0.5))
stopifnot(inherits(metage_stamp_plot(stamp), "ggplot"))
importance <- data.frame(Feature = paste0("F", 1:5), MeanDecreaseGini = 5:1)
stopifnot(inherits(metage_random_forest_importance(importance), "ggplot"))
sig <- data.frame(Feature = paste0("F", 1:5), Padj = c(0.01, 0.02, 0.04, 0.2, 0.5))
stopifnot(inherits(metage_significant_heatmap(x, sig), "ggplot"))

blob <- data.frame(GC = runif(20), Coverage = rexp(20), Bin = rep(c("B1", "B2"), 10))
stopifnot(inherits(metage_gc_coverage(blob), "ggplot"))
stopifnot(inherits(metage_bin_heatmap(x), "ggplot"))
stopifnot(inherits(metage_overlap_plot(list(A = 1:5, B = 4:8)), "ggplot"))
if (requireNamespace("ape", quietly = TRUE) && requireNamespace("ggtree", quietly = TRUE)) {
  tree_device <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tree_device)
  tree <- ape::rtree(6, tip.label = paste0("F", 1:6))
  tree_plot <- metage_tree_abundance(tree, setNames(6:1, paste0("F", 1:6)))
  stopifnot(inherits(tree_plot, "metage_tree_abundance"))
  if (requireNamespace("gridExtra", quietly = TRUE))
    stopifnot(inherits(metage_tree_abundance_grob(tree_plot), "gtable"))
  grDevices::dev.off()
  unlink(tree_device)
}

tmp <- tempfile(fileext = ".png")
metage_save_plot(metage_abundance_bar(x), tmp, width = 5, height = 4, dpi = 100)
stopifnot(file.exists(tmp), file.info(tmp)$size > 0)
unlink(tmp)
message("metage smoke tests passed")

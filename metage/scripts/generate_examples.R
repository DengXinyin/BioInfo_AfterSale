args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else file.path(tempdir(), "metage_examples")

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", all_args[grep("^--file=", all_args)])
root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
invisible(lapply(list.files(file.path(root, "metage", "R"), pattern = "[.]R$", full.names = TRUE), source))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(287)
abundance <- matrix(rpois(30 * 9, lambda = 30), 30, 9,
                    dimnames = list(paste0("Feature_", 1:30), paste0("Sample_", 1:9)))
group <- setNames(rep(c("Control", "Treatment", "Recovery"), each = 3), colnames(abundance))

base <- data.frame(Position = 1:100, A = runif(100, 22, 28), T = runif(100, 22, 28),
                   G = runif(100, 22, 28), C = runif(100, 22, 28))
metage_save_plot(metage_base_content_plot(base), file.path(output_dir, "01_base_content.png"))
metage_save_plot(metage_abundance_bar(abundance), file.path(output_dir, "02_composition.png"), 9, 6)
metage_save_plot(metage_abundance_heatmap(abundance), file.path(output_dir, "03_heatmap.png"), 9, 7)

alpha <- metage_alpha_indices(abundance)
alpha$Group <- group[alpha$Sample]
metage_save_plot(metage_alpha_boxplot(alpha, "Shannon", "Group")$plot,
                  file.path(output_dir, "04_alpha.png"))
ord <- metage_ordination(abundance, group, "PCoA")
metage_save_plot(ord$plot, file.path(output_dir, "05_pcoa.png"))
metage_save_plot(metage_diff_boxplot(abundance, group, rownames(abundance)[1], method = "anova"),
                  file.path(output_dir, "06_difference.png"))

blob <- data.frame(GC = runif(300, 0.25, 0.75), Coverage = rexp(300, 0.08),
                   Bin = sample(paste0("Bin_", 1:5), 300, TRUE))
metage_save_plot(metage_gc_coverage(blob), file.path(output_dir, "07_gc_coverage.png"))

message("Generated examples in: ", normalizePath(output_dir))

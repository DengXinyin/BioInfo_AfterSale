#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "micro/scripts/generate_examples.R"
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
micro_dir <- file.path(root, "micro")

source(file.path(micro_dir, "R", "utils.R"))
for (f in sort(list.files(file.path(micro_dir, "R"), pattern = "^[0-9].*[.]R$", full.names = TRUE))) source(f)
dir.create(file.path(micro_dir, "images"), recursive = TRUE, showWarnings = FALSE)

set.seed(542)
features <- paste0("Taxon_", sprintf("%02d", 1:35))
samples <- paste0("Sample_", sprintf("%02d", 1:12))
group <- stats::setNames(rep(c("Control", "Treatment"), each = 6), samples)
abundance <- matrix(stats::rpois(length(features) * length(samples), lambda = 50),
                    nrow = length(features), dimnames = list(features, samples))
abundance[1:6, group == "Treatment"] <- abundance[1:6, group == "Treatment"] + 100
abundance[sample(length(abundance), 100)] <- 0

p <- micro_abundance_bar(abundance, group, top_n = 12, title = "Synthetic community composition")
micro_save_plot(p, file.path(micro_dir, "images", "01_composition.png"), 9, 6)

alpha <- data.frame(Sample = samples, Group = unname(group),
                    Shannon = vegan::diversity(t(abundance), index = "shannon"))
p <- micro_alpha_boxplot(alpha, "Shannon", "Group", "Sample", title = "Shannon diversity")
micro_save_plot(p, file.path(micro_dir, "images", "02_alpha_diversity.png"), 7, 6)

ord <- micro_ordination(abundance, group, method = "PCoA", distance = "bray",
                        labels = TRUE, ellipse = TRUE)
micro_save_plot(ord$plot, file.path(micro_dir, "images", "03_ordination.png"), 8, 6)

stamp <- micro_stamp_plot(abundance, group, top_n = 12)
micro_save_plot(stamp$plot, file.path(micro_dir, "images", "04_differential.png"), 11, 7)

d <- vegan::vegdist(t(abundance), method = "bray")
anosim <- micro_anosim_plot(d, group, permutations = 199)
micro_save_plot(anosim$plot, file.path(micro_dir, "images", "05_beta_statistics.png"), 7, 6)

sets <- list(
  Control = paste0("Taxon_", sample(1:50, 28)),
  Treatment = paste0("Taxon_", sample(1:50, 31)),
  Core = paste0("Taxon_", sample(1:50, 24))
)
p <- micro_overlap_plot(sets)
micro_save_plot(p, file.path(micro_dir, "images", "06_overlap_network.png"), 7, 6)

ncm <- micro_ncm_fit(abundance)
p <- micro_ncm_plot(ncm)
micro_save_plot(p, file.path(micro_dir, "images", "07_special_plots.png"), 8, 6)

message("Generated 7 PNG files in ", file.path(micro_dir, "images"))

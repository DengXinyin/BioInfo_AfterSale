#!/usr/bin/env Rscript

# WGCNA example using the previously validated 935-gene x 15-sample dataset.
#
# Required R packages:
#   - WGCNA (direct non-base dependency)
#   - dynamicTreeCut and fastcluster (installed as WGCNA dependencies)
# Base packages used by the WGCNA_ functions:
#   - stats, graphics, grDevices, utils, tools
#
# Usage:
#   Rscript inst/examples/WGCNA_example_tutorial.R [project_dir] [output_dir]
#
# Optional environment variables:
#   BIOINFO_AFTERSALE_SOURCE=/path/to/BioInfo_AfterSale
#   WGCNA_THREADS=8
#   WGCNA_RUN_SENSITIVITY=true

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1]] else {
  "/home/xydeng/SNKH042726062901_余杰辉_转录生信分析售后"
}
output_dir <- if (length(args) >= 2L) args[[2]] else {
  file.path(tempdir(), "BioInfoAfterSale_WGCNA_Tutorial")
}
package_source <- Sys.getenv(
  "BIOINFO_AFTERSALE_SOURCE",
  unset = "/home/xydeng/Useful_Docker_Images/BioInfo_AfterSale"
)
threads <- suppressWarnings(as.integer(Sys.getenv("WGCNA_THREADS", unset = "8")))
if (is.na(threads) || threads < 1L) stop("WGCNA_THREADS must be a positive integer")

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("Package 'WGCNA' is required. Install it when building the R environment.")
}

# Load only the WGCNA module when the complete package is not installed in the
# reference image. In a complete package environment, library(BioInfoAfterSale)
# can be used instead.
source(file.path(package_source, "R", "WGCNA.R"))

expression_file <- file.path(project_dir, "Input", "selected_gene_expression_fpkm.tsv")
sample_file <- file.path(project_dir, "Input", "sampleinfo.tsv")
if (!file.exists(expression_file)) stop("Missing expression input: ", expression_file)
if (!file.exists(sample_file)) stop("Missing sample input: ", sample_file)

expression <- read.delim(
  expression_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)
sample_info <- read.delim(
  sample_file, check.names = FALSE, stringsAsFactors = FALSE
)
if (!setequal(colnames(expression), sample_info$Sample_ID)) {
  stop("Expression columns and sample_info$Sample_ID do not match")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Exact v1.6.7 image mode for this customer request: all supplied genes enter
# QC and module construction; no external phenotype is invented.
result <- WGCNA_run(
  expression,
  sample_info = sample_info,
  time_col = "time_group",
  prepare_args = list(
    transform = "none",
    keep_all_genes = TRUE,
    remove_bad = FALSE
  ),
  power_args = list(
    network_type = "unsigned",
    strategy = "image_rule",
    fit_cutoff = 0.9
  ),
  module_args = list(
    deep_split = 2,
    min_module_size = 30,
    merge_cut_height = 0.25,
    reassign_threshold = 0,
    pam_respects_dendro = FALSE,
    max_block_size = 1000,
    threads = threads,
    seed = 123
  ),
  output_dir = output_dir
)

WGCNA_plot_soft_threshold(
  result$power,
  file.path(output_dir, "Soft_thresholding_power.pdf")
)

# Export the strongest network edges from the largest non-grey module.
summary_table <- result$modules$module_summary
eligible <- summary_table[summary_table$module != "grey", , drop = FALSE]
largest_module <- eligible$module[[which.max(eligible$gene_count)]]
network_result <- WGCNA_export_network(
  result$modules,
  module = largest_module,
  threshold = 0.15,
  max_module_genes = 2000,
  max_edges = 5000,
  output_dir = file.path(output_dir, "Network")
)

# Regression check for the exact previously validated dataset.
if (nrow(expression) == 935L && ncol(expression) == 15L) {
  expected <- c(
    turquoise = 316, blue = 238, brown = 132, yellow = 87,
    green = 73, red = 72, grey = 17
  )
  observed <- setNames(summary_table$gene_count, summary_table$module)
  if (!isTRUE(all.equal(
    as.numeric(observed[names(expected)]), as.numeric(expected)
  ))) {
    warning("Module sizes differ from the validated raw-FPKM/image-rule result")
  }
}

# Optional scale sensitivity analysis. This is intentionally opt-in because it
# performs a second complete network construction.
if (tolower(Sys.getenv("WGCNA_RUN_SENSITIVITY", unset = "false")) == "true") {
  WGCNA_run(
    expression,
    sample_info = sample_info,
    time_col = "time_group",
    prepare_args = list(
      transform = "log2p1",
      keep_all_genes = TRUE,
      remove_bad = FALSE
    ),
    power_args = list(
      network_type = "unsigned",
      strategy = "scale_free",
      fit_cutoff = 0.9
    ),
    module_args = list(threads = threads, seed = 123),
    output_dir = file.path(output_dir, "log2p1_sensitivity"),
    verbose = TRUE
  )
}

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
saveRDS(result, file.path(output_dir, "WGCNA_workflow_result.rds"))

cat("WGCNA tutorial completed\n")
cat("  Project:", normalizePath(project_dir), "\n")
cat("  Output:", normalizePath(output_dir), "\n")
cat("  Power:", result$power$power, "\n")
cat("  Modules including grey:", nrow(summary_table), "\n")
cat("  Genes assigned:", sum(summary_table$gene_count), "\n")

\name{WGCNA_workflow}
\alias{WGCNA_prepare_expression}
\alias{WGCNA_select_power}
\alias{WGCNA_build_modules}
\alias{WGCNA_module_sample}
\alias{WGCNA_module_trait}
\alias{WGCNA_module_membership}
\alias{WGCNA_export_network}
\alias{WGCNA_run}
\alias{WGCNA_write_results}
\alias{WGCNA_plot_soft_threshold}
\title{Reusable WGCNA workflow functions}
\description{
Modular functions for expression quality control, soft-threshold selection,
module construction, module-sample or module-trait summaries, gene-module
membership, TOM network export, result writing, and diagnostic plotting.
Defaults are derived from the v1.6.7 WGCNA workflow image, while key
statistical choices remain explicit and configurable.
}
\usage{
WGCNA_prepare_expression(expression, gene_id_col = NULL,
  samples_in_rows = FALSE, transform = c("none", "log2p1"),
  zero_fraction_max = 0.3, keep_all_genes = FALSE,
  remove_bad = TRUE, verbose = 0)

WGCNA_select_power(prepared,
  powers = c(1:10, seq(12, 20, by = 2)),
  network_type = c("unsigned", "signed", "signed hybrid"),
  strategy = c("scale_free", "image_rule", "manual"),
  fit_cutoff = 0.9, manual_power = NULL, verbose = 0)

WGCNA_build_modules(prepared, power, network_type = NULL,
  TOM_type = NULL, deep_split = 2, min_module_size = 30,
  merge_cut_height = 0.25, reassign_threshold = 0,
  pam_respects_dendro = FALSE, max_block_size = 1000,
  save_TOMs = FALSE, threads = 1, seed = 123, verbose = 0)

WGCNA_module_sample(modules, sample_info = NULL, time_col = NULL,
  orient = TRUE)

WGCNA_module_trait(modules, traits, MEs = NULL,
  use = "pairwise.complete.obs", p_adjust_method = "bonferroni")

WGCNA_module_membership(modules, MEs = NULL,
  p_adjust_method = "bonferroni")

WGCNA_export_network(modules, module = NULL, threshold = 0.15,
  max_module_genes = 2000, max_edges = Inf, include_grey = FALSE,
  output_dir = NULL)

WGCNA_run(expression, sample_info = NULL, traits = NULL,
  prepare_args = list(), power_args = list(), module_args = list(),
  time_col = NULL, output_dir = NULL, export_networks = FALSE,
  network_args = list(), verbose = TRUE)

WGCNA_write_results(result, output_dir)

WGCNA_plot_soft_threshold(power_result, output_file = NULL,
  width = 16, height = 8)
}
\arguments{
  \item{expression}{Numeric expression matrix, normally genes by samples.}
  \item{gene_id_col}{Optional expression column containing gene IDs.}
  \item{samples_in_rows}{Whether samples are already rows.}
  \item{transform}{Optional \code{log2p1} transformation.}
  \item{zero_fraction_max}{Maximum zero fraction retained per gene.}
  \item{keep_all_genes}{Disable the zero-fraction filter.}
  \item{remove_bad}{Remove or reject data flagged by \code{goodSamplesGenes}.}
  \item{verbose}{WGCNA or workflow verbosity.}
  \item{prepared}{Prepared object or sample-by-gene matrix.}
  \item{powers}{Candidate soft powers.}
  \item{network_type}{Unsigned, signed, or signed-hybrid network.}
  \item{strategy}{Soft-power selection strategy.}
  \item{fit_cutoff}{Scale-free topology cutoff.}
  \item{manual_power}{Explicit power for the manual strategy.}
  \item{power}{Numeric power or a power-selection object.}
  \item{TOM_type}{TOM type; defaults to the network type.}
  \item{deep_split,min_module_size,merge_cut_height,reassign_threshold,
  pam_respects_dendro,max_block_size,save_TOMs,threads,seed}{Module-building parameters.}
  \item{modules}{Object returned by \code{WGCNA_build_modules}.}
  \item{sample_info}{Sample metadata with aligned IDs.}
  \item{time_col}{Metadata column used for descriptive grouped summaries.}
  \item{orient}{Orient eigengenes toward module mean expression.}
  \item{traits}{Numeric sample-by-trait data frame.}
  \item{MEs}{Optional eigengene matrix.}
  \item{use}{Correlation missing-value policy.}
  \item{p_adjust_method}{Multiple-testing correction method.}
  \item{module}{Modules to export.}
  \item{threshold,max_module_genes,max_edges,include_grey}{Network export controls.}
  \item{output_dir}{Optional result directory.}
  \item{prepare_args,power_args,module_args,network_args}{Named argument lists.}
  \item{export_networks}{Whether to export TOM networks.}
  \item{result}{Workflow result returned by \code{WGCNA_run}.}
  \item{power_result}{Power-selection result.}
  \item{output_file}{Optional PDF or PNG diagnostic plot.}
  \item{width,height}{Plot dimensions in inches.}
}
\value{
Each step returns a named list containing both reusable intermediate objects
and tidy result tables. \code{WGCNA_run} returns all workflow components
invisibly.
}
\details{
For no-trait analyses, module-sample relationships are represented by module
eigengene values and optional grouped summaries. The workflow intentionally
does not create one-sample dummy traits or report pseudo-correlation P values.
}

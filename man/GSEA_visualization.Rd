\name{GSEA_visualization}
\alias{GSEA_visualization}
\alias{calculate_gsea_curve}
\alias{plot_gsea}
\alias{plot_gsea_pathway}
\alias{plot_gsea_result}
\alias{plot_gsea_batch}
\title{Calculate and visualize GSEA running scores}
\usage{
calculate_gsea_curve(ranked_metric, gene_set, exponent = 1,
  sort = TRUE, na_rm = TRUE)

plot_gsea(running_score, hits, ranked_metric, title = NULL,
  statistics = NULL, colors = c("#E25659", "#335372"),
  metric_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
  show_metric_band = TRUE, peak_position = NULL, style = NULL,
  output_file = NULL, width = NULL, height = NULL, dpi = NULL)

plot_gsea_pathway(ranked_metric, gene_set, exponent = 1,
  title = NULL, statistics = NULL,
  colors = c("#E25659", "#335372"),
  metric_colors = c("#2166AC", "#F7F7F7", "#B2182B"),
  show_metric_band = TRUE, style = NULL, output_file = NULL,
  width = NULL, height = NULL, dpi = NULL)

plot_gsea_result(result, gene_set_id, title = NULL, ...)

plot_gsea_batch(result, gene_set_ids, output_dir,
  format = c("pdf", "png"), filename_prefix = "GSEA",
  width = 9.2, height = 6.8, dpi = 300, ...)
}
\arguments{
  \item{ranked_metric}{Named numeric ranked gene statistic.}
  \item{gene_set}{Character vector of genes in one pathway.}
  \item{exponent}{Non-negative GSEA weighting exponent.}
  \item{sort}{Whether to sort the metric decreasingly.}
  \item{na_rm}{Whether to remove non-finite metric values.}
  \item{running_score}{Numeric running enrichment score.}
  \item{hits}{Integer hit positions or a logical vector.}
  \item{title}{Optional pathway title.}
  \item{statistics}{Optional named statistics shown in the top panel.}
  \item{colors}{Positive and negative plot colors.}
  \item{metric_colors}{Three colors for the ranked-metric band.}
  \item{show_metric_band}{Whether to draw the middle color band.}
  \item{peak_position}{Optional rank of the ES peak.}
  \item{style}{A \code{bioinfo_plot_style} object.}
  \item{output_file}{Optional PDF or PNG output path.}
  \item{width, height, dpi}{Output dimensions and resolution.}
  \item{result}{A \code{clusterProfiler::gseaResult} object.}
  \item{gene_set_id}{One exact pathway ID or Description.}
  \item{gene_set_ids}{Pathway IDs or Descriptions to export.}
  \item{output_dir}{Batch output directory.}
  \item{format}{Batch output format.}
  \item{filename_prefix}{Batch filename prefix.}
  \item{...}{Additional plotting arguments.}
}
\value{
\code{calculate_gsea_curve()} returns a \code{bioinfo_gsea_curve} list.
The plotting functions return ggplot objects. \code{plot_gsea_batch()}
invisibly returns the generated file paths.
}
\description{
Calculates weighted GSEA running scores and draws a three-panel pathway figure
with the enrichment curve, hit positions, ranked-metric color band, and ranked
metric distribution. Convenience wrappers support clusterProfiler GSEA objects
and batch export.
}
\examples{
set.seed(1)
metric <- sort(setNames(rnorm(100), paste0("G", 1:100)), decreasing = TRUE)
genes <- names(metric)[1:15]
curve <- calculate_gsea_curve(metric, genes)
}

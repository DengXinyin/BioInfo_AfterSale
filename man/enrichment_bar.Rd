\name{enrichment_bar}
\alias{enrichment_bar}
\title{Advanced Enrichment Bar Plot}
\usage{
enrichment_bar(
  data,
  style,
  output_file = NULL,
  width = NULL,
  height = NULL,
  x_axis_label = "pvalue",
  component_spacing = 0.4,
  group_label_scale = 1,
  bubble_label_scale = 0.6,
  description_scale = 1,
  gene_scale = 1
)
}
\arguments{
\item{data}{A data.frame, or a path to an Excel file. Required columns are
\code{group}, \code{Description}, \code{pvalue}, \code{count}, and
\code{geneID}.}

\item{style}{Object returned by \code{\link{choose_plot_style}}. Its global
settings, text styles, legend options, panel theme, and group palette control
the plot appearance.}

\item{output_file}{Optional PNG output path. If \code{NULL}, the ggplot object
is returned without writing a file.}

\item{width, height}{Output size in inches. \code{NULL} uses
\code{style$global$figure_width} and \code{style$global$figure_height}.}

\item{x_axis_label}{Text displayed inside \code{-log10()} on the x axis.}

\item{component_spacing}{Distance between the group strip, count bubble, and
bar.}

\item{group_label_scale, bubble_label_scale, description_scale, gene_scale}{
Multipliers applied after the text sizes supplied by \code{style}. Set
\code{gene_scale = 0} to hide gene labels.}
}
\value{
A ggplot object, invisibly when \code{output_file} is supplied.
}
\description{
Draws a grouped enrichment bar plot whose visual style is managed by
\code{\link{choose_plot_style}}.
}

#' Plot a three-panel GSEA running-score figure
#'
#' @param running_score Numeric running enrichment score.
#' @param hits Integer hit positions or a logical vector.
#' @param ranked_metric Numeric ranked metric of the same length.
#' @param title Optional title.
#' @param statistics Optional named values displayed in the top panel.
#' @param colors Positive and negative colors.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A faceted ggplot object.
#' @export
plot_gsea <- function(running_score, hits, ranked_metric, title = NULL,
                      statistics = NULL, colors = c("#E25659", "#335372"),
                      style = NULL, output_file = NULL, width = NULL,
                      height = NULL, dpi = NULL) {
  style <- .plot_style_or_default(style); n <- length(running_score)
  if (length(ranked_metric) != n) stop("`running_score` and `ranked_metric` must have equal lengths.", call. = FALSE)
  if (is.logical(hits)) hits <- which(hits)
  if (any(hits < 1 | hits > n)) stop("`hits` contains out-of-range positions.", call. = FALSE)
  panels <- factor(c(rep("Running Enrichment Score",n),rep("Hits",length(hits)),rep("Ranked Metric",n)), levels=c("Running Enrichment Score","Hits","Ranked Metric"))
  d <- data.frame(x=c(seq_len(n),hits,seq_len(n)), y=c(running_score,rep(.5,length(hits)),ranked_metric), panel=panels)
  p <- ggplot2::ggplot(d,ggplot2::aes(.data$x,.data$y)) +
    ggplot2::geom_line(data=d[d$panel=="Running Enrichment Score",],color=colors[1],linewidth=.8) +
    ggplot2::geom_segment(data=d[d$panel=="Hits",],ggplot2::aes(xend=.data$x,y=0,yend=1),color=colors[2],linewidth=.25) +
    ggplot2::geom_area(data=d[d$panel=="Ranked Metric" & d$y>=0,],fill=colors[1],alpha=.75) +
    ggplot2::geom_area(data=d[d$panel=="Ranked Metric" & d$y<0,],fill=colors[2],alpha=.75) +
    ggplot2::facet_grid(rows=ggplot2::vars(.data$panel),scales="free_y",switch="y") +
    ggplot2::labs(x="Rank in Ordered Gene List",y=NULL,title=title) + style$ggplot_theme +
    ggplot2::theme(panel.spacing.y=grid::unit(0,"pt"),strip.placement="outside",strip.background=ggplot2::element_blank())
  if (!is.null(statistics)) {
    label <- paste(names(statistics),format(statistics,digits=3),sep=" = ",collapse="\n")
    p <- p + ggplot2::annotate("text",x=Inf,y=Inf,label=label,hjust=1.05,vjust=1.15)
  }
  .save_bioinfo_plot(p, style, output_file, width, height, dpi); p
}

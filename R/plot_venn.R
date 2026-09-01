#' Plot a three-set Venn diagram
#'
#' @param sets Named list of exactly three vectors.
#' @param colors Optional three colors; defaults to the style palette.
#' @param title Optional title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A ggplot object.
#' @export
plot_venn <- function(sets, colors = NULL, title = NULL, style = NULL,
                      output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  style <- .plot_style_or_default(style)
  if (!is.list(sets) || length(sets) != 3L || is.null(names(sets))) stop("`sets` must be a named list of three vectors.", call. = FALSE)
  sets <- lapply(sets, unique); if (is.null(colors)) colors <- style$group_palette[1:3]
  circle <- function(cx, cy, r = 1, n = 300) data.frame(x = cx + r*cos(seq(0, 2*pi, length.out=n)), y = cy + r*sin(seq(0, 2*pi, length.out=n)))
  xy <- matrix(c(-.7,.35,.7,.35,0,-.55), ncol=2, byrow=TRUE)
  poly <- do.call(rbind, lapply(1:3, function(i) transform(circle(xy[i,1],xy[i,2]), set=names(sets)[i])))
  a <- sets[[1]]; b <- sets[[2]]; c <- sets[[3]]
  counts <- c(length(setdiff(a,union(b,c))), length(setdiff(b,union(a,c))), length(setdiff(c,union(a,b))),
              length(setdiff(intersect(a,b),c)), length(setdiff(intersect(a,c),b)), length(setdiff(intersect(b,c),a)), length(Reduce(intersect,sets)))
  pos <- data.frame(x=c(-1.05,1.05,0,-.0,-.55,.55,0), y=c(.45,.45,-.95,.65,-.15,-.15,.05), label=counts)
  pal <- stats::setNames(colors, names(sets))
  p <- ggplot2::ggplot(poly, ggplot2::aes(.data$x,.data$y,group=.data$set,fill=.data$set)) +
    ggplot2::geom_polygon(alpha=.42,color=NA) + ggplot2::geom_text(data=pos,ggplot2::aes(.data$x,.data$y,label=.data$label),inherit.aes=FALSE) +
    ggplot2::annotate("text",x=xy[,1],y=c(1.5,1.5,-1.65),label=names(sets)) +
    ggplot2::scale_fill_manual(values=pal) + ggplot2::coord_equal() + ggplot2::labs(title=title,fill=NULL) +
    style$ggplot_theme + ggplot2::theme(axis.title=ggplot2::element_blank(),axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),panel.grid=ggplot2::element_blank())
  .save_bioinfo_plot(p, style, output_file, width, height, dpi); p
}

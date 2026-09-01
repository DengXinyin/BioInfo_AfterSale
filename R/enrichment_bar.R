#!/usr/bin/env Rscript

# Advanced enrichment bar plot compatible with BioInfo_AfterSale::choose_plot_style().
# Required columns: group, Description, pvalue, count, geneID
# All enrichment-bar text prefers Times New Roman through the supplied style.

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
})

.bar_text_face <- function(style) {
  if (isTRUE(style$bold) && isTRUE(style$italic)) return("bold.italic")
  if (isTRUE(style$bold)) return("bold")
  if (isTRUE(style$italic)) return("italic")
  "plain"
}

.bar_text_family <- function(style, global_family) {
  if (!is.null(style$font_family) && nzchar(style$font_family)) {
    style$font_family
  } else {
    global_family
  }
}

.bar_text_element <- function(style, global_family) {
  if (!isTRUE(style$show)) return(ggplot2::element_blank())
  ggplot2::element_text(
    family = .bar_text_family(style, global_family),
    size = style$size,
    face = .bar_text_face(style),
    hjust = c(left = 0, center = 0.5, right = 1)[[style$align]],
    color = "black"
  )
}

.bar_portable_family <- function(font_family) {
  if (tolower(font_family) == "arial") {
    "sans"
  } else if (tolower(font_family) == "times new roman") {
    "serif"
  } else {
    font_family
  }
}

#' Draw an advanced enrichment bar plot with a shared BioInfo_AfterSale style.
#'
#' @param data A data.frame, or a path to an Excel file. Required columns:
#'   `group`, `Description`, `pvalue`, `count`, and `geneID`.
#' @param style Object returned by `choose_plot_style()` from the
#'   BioInfo_AfterSale repository. Load that function from the repository before
#'   calling this function.
#' @param output_file Optional PNG output path. If NULL, the ggplot object is returned.
#' @param width,height Output size in inches. NULL uses `style$global`.
#' @param x_axis_label Text displayed inside `-log10()`.
#' @param component_spacing Distance between the group strip, count bubble and bar.
#' @param group_label_scale,bubble_label_scale,description_scale,gene_scale
#'   Multipliers applied after the text sizes supplied by `style`.
#' @return A ggplot object invisibly when `output_file` is supplied.
enrichment_bar <- function(
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
    gene_scale = 1) {

  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
  }
  if (is.null(width)) width <- style$global$figure_width
  if (is.null(height)) height <- style$global$figure_height
  if (component_spacing <= 0) stop("component_spacing must be greater than 0.")

  if (is.character(data) && length(data) == 1L) {
    data <- as.data.frame(suppressMessages(suppressWarnings(read_excel(data))))
  }
  if (!is.data.frame(data)) stop("data must be a data.frame or an Excel-file path.")

  required_columns <- c("group", "Description", "pvalue", "count", "geneID")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  dat <- data %>%
    transmute(
      group = as.character(group),
      Description = as.character(Description),
      pvalue = as.numeric(pvalue),
      count = as.numeric(count),
      geneID = gsub("/", "; ", as.character(geneID))
    ) %>%
    filter(!is.na(group), !is.na(Description), is.finite(pvalue), pvalue > 0,
           is.finite(count))
  if (nrow(dat) == 0) stop("No valid enrichment records were found.")

  dat$group <- factor(dat$group, levels = unique(dat$group))
  dat$score <- -log10(dat$pvalue)
  dat <- dat %>%
    mutate(row = row_number(),
           y = rev(row) + (nlevels(group) - as.integer(group)) * 0.4)
  group_ranges <- dat %>%
    group_by(group) %>%
    summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, ymid = mean(y), .groups = "drop")

  group_colours <- setNames(
    rep(style$group_palette, length.out = nlevels(dat$group)),
    levels(dat$group)
  )

  x_bubble <- -component_spacing
  x_strip_right <- -(component_spacing + 0.4)
  x_left <- x_strip_right - 0.65
  x_max <- max(dat$score) * 1.06
  gene_y <- dat$y - 0.48

  global_family <- .bar_portable_family(style$global$font_family)
  data_label <- style$text$data_label
  description_label <- style$text$axis_text
  gene_label <- style$text$data_label
  group_label <- style$text$facet_label

  plot <- ggplot(dat, aes(y = y)) +
    geom_rect(data = group_ranges,
      aes(xmin = x_left, xmax = x_strip_right, ymin = ymin, ymax = ymax, fill = group),
      inherit.aes = FALSE, alpha = 0.33, colour = NA, show.legend = FALSE) +
    geom_rect(aes(xmin = 0, xmax = score, ymin = y - 0.30, ymax = y + 0.30, fill = group),
      colour = NA, show.legend = FALSE) +
    geom_point(aes(x = x_bubble, size = count, fill = group), shape = 21,
      colour = "black", stroke = 0.75, show.legend = TRUE) +
    geom_text(aes(x = x_bubble, label = format(count, trim = TRUE)),
      size = data_label$size / 3.2 * bubble_label_scale,
      family = .bar_text_family(data_label, global_family),
      fontface = .bar_text_face(data_label)) +
    geom_text(aes(x = 0.10, label = Description), hjust = 0,
      size = description_label$size / 3.2 * description_scale,
      family = .bar_text_family(description_label, global_family),
      fontface = .bar_text_face(description_label), colour = "black") +
    geom_text(aes(x = 0.10, y = gene_y, label = geneID, colour = group), hjust = 0,
      size = gene_label$size / 3.2 * gene_scale,
      family = .bar_text_family(gene_label, global_family),
      fontface = .bar_text_face(gene_label), show.legend = FALSE) +
    geom_text(data = group_ranges,
      aes(x = mean(c(x_left, x_strip_right)), y = ymid, label = group),
      inherit.aes = FALSE, angle = 90,
      size = group_label$size / 3.2 * group_label_scale,
      family = .bar_text_family(group_label, global_family),
      fontface = .bar_text_face(group_label)) +
    scale_fill_manual(values = group_colours, name = "Group") +
    scale_colour_manual(values = group_colours, guide = "none") +
    scale_size_continuous(name = "Gene Count", range = c(5, 15),
      breaks = unique(round(seq(min(dat$count), max(dat$count), length.out = 3)))) +
    scale_x_continuous(limits = c(x_left, x_max), breaks = seq(0, floor(x_max), by = 2),
      expand = expansion(mult = c(0, 0.06))) +
    scale_y_continuous(limits = c(0.15, max(dat$y) + 0.55), breaks = NULL, expand = c(0, 0)) +
    labs(x = bquote(-log[10](.(x_axis_label))), y = NULL) +
    guides(size = guide_legend(order = 1, override.aes = list(fill = "white", colour = "black")),
      fill = guide_legend(order = 2, override.aes = list(shape = 22, size = 5))) +
    coord_cartesian(clip = "off") +
    style$ggplot_theme +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = .bar_text_element(style$text$axis_title, global_family),
      axis.text.x = .bar_text_element(style$text$axis_text, global_family),
      legend.title = .bar_text_element(style$text$legend_title, global_family),
      legend.text = .bar_text_element(style$text$legend_text, global_family),
      legend.box.spacing = unit(0.2, "cm"),
      plot.margin = margin(16, 10, 24, 40)
    )

  remove_layers <- integer()
  if (!isTRUE(style$text$data_label$show)) remove_layers <- c(remove_layers, 4L)
  if (gene_scale == 0) remove_layers <- c(remove_layers, 6L)
  if (!isTRUE(style$text$facet_label$show)) remove_layers <- c(remove_layers, 7L)
  if (length(remove_layers)) {
    plot$layers <- plot$layers[-sort(unique(remove_layers), decreasing = TRUE)]
  }

  if (!is.null(output_file)) {
    ggsave(output_file, plot, width = width, height = height, dpi = style$global$dpi, bg = "white")
    message("Wrote: ", normalizePath(output_file))
  }
  invisible(plot)
}

# All plot text prefers Times New Roman; this is also the default family.
.resolve_plot_font_family <- function(font_family) {
  if (!identical(font_family, "Times New Roman") || !nzchar(Sys.which("fc-match"))) {
    return(font_family)
  }
  matched <- tryCatch(
    paste(system2("fc-match", c("-f", "%{family}", font_family), stdout = TRUE), collapse = " "),
    error = function(e) ""
  )
  if (grepl("Times New Roman", matched, fixed = TRUE)) font_family else "serif"
}

.plot_text_defaults <- function() {
  list(
    title = list(font_family = "", size = 24, bold = TRUE, italic = FALSE,
                 align = "center", show = TRUE),
    subtitle = list(font_family = "", size = 20, bold = FALSE, italic = FALSE,
                    align = "center", show = TRUE),
    axis_title = list(font_family = "", size = 20, bold = FALSE, italic = FALSE,
                      align = "center", show = TRUE),
    axis_text = list(font_family = "", size = 18, bold = FALSE, italic = FALSE,
                     align = "center", show = TRUE),
    legend_title = list(font_family = "", size = 20, bold = FALSE, italic = FALSE,
                        align = "left", show = TRUE),
    legend_text = list(font_family = "", size = 18, bold = FALSE, italic = FALSE,
                       align = "left", show = TRUE),
    data_label = list(font_family = "", size = 16, bold = FALSE, italic = FALSE,
                      align = "center", show = TRUE),
    facet_label = list(font_family = "", size = 18, bold = TRUE, italic = FALSE,
                       align = "center", show = TRUE)
  )
}

.plot_group_palette <- function() {
  c(
    "#4472C4", "#ED7D31", "#70AD47", "#A5A5A5", "#FFC000", "#5B9BD5",
    "#C55A11", "#8064A2", "#2F5597", "#A9D18E", "#F4B183", "#9E480E"
  )
}

.merge_plot_style <- function(base, override) {
  if (!is.list(override)) stop("Plot-style overrides must be lists.", call. = FALSE)
  if (!length(override)) return(base)
  if (is.null(names(override)) || any(!nzchar(names(override)))) {
    stop("Plot-style overrides must be named lists.", call. = FALSE)
  }
  unknown <- setdiff(names(override), names(base))
  if (length(unknown)) {
    stop("Unknown plot-style field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  for (name in names(override)) base[[name]] <- override[[name]]
  base
}

.validate_text_style <- function(style, element) {
  if (!is.character(style$font_family) || length(style$font_family) != 1L ||
      is.na(style$font_family)) {
    stop("`", element, "$font_family` must be one character value.", call. = FALSE)
  }
  .assert_positive_number(style$size, paste0(element, "$size"))
  for (field in c("bold", "italic", "show")) {
    if (!is.logical(style[[field]]) || length(style[[field]]) != 1L ||
        is.na(style[[field]])) {
      stop("`", element, "$", field, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  .assert_choice(style$align, c("left", "center", "right"), paste0(element, "$align"))
  style
}

.text_face <- function(style) {
  if (style$bold && style$italic) return("bold.italic")
  if (style$bold) return("bold")
  if (style$italic) return("italic")
  "plain"
}

.text_element <- function(style, global_family) {
  if (!style$show) return(ggplot2::element_blank())
  family <- if (nzchar(style$font_family)) style$font_family else global_family
  ggplot2::element_text(
    family = family,
    size = style$size,
    face = .text_face(style),
    hjust = c(left = 0, center = 0.5, right = 1)[[style$align]],
    color = "black"
  )
}

.build_ggplot_theme <- function(config) {
  base <- switch(
    config$global$theme,
    bw = ggplot2::theme_bw(base_family = config$global$font_family),
    classic = ggplot2::theme_classic(base_family = config$global$font_family)
  )
  text <- config$text
  legend <- config$legend
  panel <- config$panel
  legend_position <- if (legend$show) legend$position else "none"
  legend_background <- if (legend$frame) {
    ggplot2::element_rect(color = "black", fill = "white", linewidth = 0.5)
  } else {
    ggplot2::element_blank()
  }

  base + ggplot2::theme(
    text = ggplot2::element_text(family = config$global$font_family),
    plot.title = .text_element(text$title, config$global$font_family),
    plot.subtitle = .text_element(text$subtitle, config$global$font_family),
    axis.title = .text_element(text$axis_title, config$global$font_family),
    axis.text = .text_element(text$axis_text, config$global$font_family),
    legend.title = .text_element(text$legend_title, config$global$font_family),
    legend.text = .text_element(text$legend_text, config$global$font_family),
    strip.text = .text_element(text$facet_label, config$global$font_family),
    legend.position = legend_position,
    legend.box = legend$box,
    legend.box.just = legend$box_just,
    legend.background = legend_background,
    legend.box.background = legend_background,
    panel.border = if (panel$border) {
      ggplot2::element_rect(
        color = panel$border_color, fill = NA, linewidth = panel$border_width
      )
    } else {
      ggplot2::element_blank()
    },
    panel.grid.major = if (panel$major_grid) {
      ggplot2::element_line(
        color = panel$major_grid_color, linewidth = panel$major_grid_width
      )
    } else {
      ggplot2::element_blank()
    },
    panel.grid.minor = if (panel$minor_grid) {
      ggplot2::element_line(
        color = panel$minor_grid_color, linewidth = panel$minor_grid_width
      )
    } else {
      ggplot2::element_blank()
    }
  )
}

#' Choose a shared ggplot2 style
#'
#' Creates a reusable visualization-style object modelled on the common style
#' interface used by the metagenomics workflow. It separates global settings,
#' semantic text elements, legend settings, and a group-color palette. The
#' object can be passed to [GO_KEGG_plot()] and future plotting functions.
#'
#' @param font_family Global font family.
#' @param theme Base ggplot2 theme: `"bw"` or `"classic"`.
#' @param dpi Output resolution stored with the style object.
#' @param figure_width Default output width in inches.
#' @param figure_height Default output height in inches.
#' @param title,subtitle,axis_title,axis_text,legend_title,legend_text,data_label,facet_label
#'   Named lists overriding `font_family`, `size`, `bold`, `italic`, `align`,
#'   or `show` for each semantic text element.
#' @param legend Named list overriding `show`, `position`, `frame`, `box`, or
#'   `box_just`.
#' @param panel Named list controlling `border`, `border_color`,
#'   `border_width`, `major_grid`, `minor_grid`, their colors, and their line
#'   widths.
#' @param group_palette Character vector of group colors. `NULL` uses the
#'   package's default 12-color palette.
#'
#' @return A `bioinfo_plot_style` list containing `global`, `text`, `legend`,
#'   `panel`, `group_palette`, and the generated `ggplot_theme`.
#' @export
choose_plot_style <- function(
    font_family = "Times New Roman",
    theme = c("bw", "classic"),
    dpi = 300,
    figure_width = 10,
    figure_height = 8,
    title = list(),
    subtitle = list(),
    axis_title = list(),
    axis_text = list(),
    legend_title = list(),
    legend_text = list(),
    data_label = list(),
    facet_label = list(),
    legend = list(),
    panel = list(),
    group_palette = NULL) {

  theme <- match.arg(theme)
  if (!is.character(font_family) || length(font_family) != 1L ||
      is.na(font_family) || !nzchar(font_family)) {
    stop("`font_family` must be one non-empty character value.", call. = FALSE)
  }
  .assert_positive_number(dpi, "dpi")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  font_family <- .resolve_plot_font_family(font_family)
  if (is.null(group_palette)) group_palette <- .plot_group_palette()

  overrides <- list(
    title = title, subtitle = subtitle, axis_title = axis_title,
    axis_text = axis_text, legend_title = legend_title,
    legend_text = legend_text, data_label = data_label,
    facet_label = facet_label
  )
  text <- .plot_text_defaults()
  for (element in names(text)) {
    text[[element]] <- .merge_plot_style(text[[element]], overrides[[element]])
    text[[element]] <- .validate_text_style(text[[element]], element)
  }

  legend <- .merge_plot_style(
    list(
      show = TRUE, position = "right", frame = FALSE,
      box = "vertical", box_just = "center"
    ),
    legend
  )
  for (field in c("show", "frame")) {
    if (!is.logical(legend[[field]]) || length(legend[[field]]) != 1L ||
        is.na(legend[[field]])) {
      stop("`legend$", field, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  .assert_choice(
    legend$position, c("left", "right", "top", "bottom", "none"),
    "legend$position"
  )
  .assert_choice(legend$box, c("vertical", "horizontal"), "legend$box")
  .assert_choice(legend$box_just, c("left", "center", "right"), "legend$box_just")

  panel <- .merge_plot_style(
    list(
      border = identical(theme, "bw"),
      border_color = "black",
      border_width = 0.8,
      major_grid = identical(theme, "bw"),
      minor_grid = FALSE,
      major_grid_color = "#D9D9D9",
      minor_grid_color = "#EEEEEE",
      major_grid_width = 0.5,
      minor_grid_width = 0.3
    ),
    panel
  )
  for (field in c("border", "major_grid", "minor_grid")) {
    .assert_flag(panel[[field]], paste0("panel$", field))
  }
  for (field in c("border_color", "major_grid_color", "minor_grid_color")) {
    value <- panel[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("`panel$", field, "` must be one non-empty color value.", call. = FALSE)
    }
  }
  for (field in c("border_width", "major_grid_width", "minor_grid_width")) {
    .assert_positive_number(panel[[field]], paste0("panel$", field))
  }
  if (!is.character(group_palette) || !length(group_palette) ||
      anyNA(group_palette) || any(!nzchar(group_palette))) {
    stop("`group_palette` must be a non-empty character vector.", call. = FALSE)
  }

  config <- list(
    global = list(
      font_family = font_family,
      theme = theme,
      dpi = dpi,
      figure_width = figure_width,
      figure_height = figure_height
    ),
    text = text,
    legend = legend,
    panel = panel,
    group_palette = group_palette
  )
  config$ggplot_theme <- .build_ggplot_theme(config)
  class(config) <- c("bioinfo_plot_style", "list")
  config
}

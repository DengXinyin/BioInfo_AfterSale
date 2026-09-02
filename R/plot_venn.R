#' Draw a three-set Venn diagram
#'
#' Builds a Venn diagram from either a data frame of 0/1 membership indicators
#' or a list of character sets, using [ggvenn::ggvenn()]. Two to three sets are
#' supported, matching the shared toolbox colour convention by default.
#'
#' When `data` is a data frame, `set_columns` names the indicator columns (each
#' value 1 means the row belongs to that set) and `id_column` names the gene or
#' feature identifier column. When `id_column` is `NULL`, row names are used.
#'
#' @param data Data frame of 0/1 membership indicators, or a named list of
#'   character vectors (one vector per set). A list input skips the data-frame
#'   indicator parsing and is passed straight to [ggvenn::ggvenn()].
#' @param set_columns Character vector naming the indicator columns. `NULL`
#'   (default) uses every column other than `id_column`.
#' @param id_column Name of the feature identifier column used to split the
#'   sets. `NULL` (default) uses row names.
#' @param set_names Optional character vector of display names, one per set in
#'   the order of `set_columns`. `NULL` uses the column names.
#' @param set_colors Named or unnamed colour vector passed to
#'   [ggvenn::ggvenn()] as `fill_color`. `NULL` uses the shared toolbox palette.
#' @param fill_alpha Fill opacity passed to [ggvenn::ggvenn()].
#' @param stroke_size Set-outline stroke width passed to
#'   [ggvenn::ggvenn()].
#' @param set_name_size Size of the set-name text.
#' @param title Optional plot title.
#' @param font_family Font family used by the shared plot style. Ignored when
#'   `style` is supplied.
#' @param style Optional object returned by [choose_plot_style()].
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object from [ggvenn::ggvenn()].
#' @export
plot_venn <- function(
    data,
    set_columns = NULL,
    id_column = NULL,
    set_names = NULL,
    set_colors = NULL,
    fill_alpha = 0.5,
    stroke_size = 0.5,
    set_name_size = 5,
    title = NULL,
    font_family = "sans",
    style = NULL,
    output_file = NULL,
    figure_width = 8,
    figure_height = 8,
    dpi = 300) {

  .assert_probability(fill_alpha, "fill_alpha")
  .assert_positive_number(stroke_size, "stroke_size")
  .assert_positive_number(set_name_size, "set_name_size")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")

  if (is.data.frame(data)) {
    if (!is.null(id_column) &&
        (!is.character(id_column) || length(id_column) != 1L || is.na(id_column))) {
      stop("`id_column` must be NULL or one column name.", call. = FALSE)
    }
    if (!is.null(id_column) && !id_column %in% names(data)) {
      stop("`id_column` is not a column in `data`.", call. = FALSE)
    }
    id <- if (is.null(id_column)) rownames(data) else as.character(data[[id_column]])
    if (is.null(set_columns)) {
      set_columns <- setdiff(names(data), id_column)
    }
    if (!is.character(set_columns) || !length(set_columns) ||
        anyNA(set_columns) || any(!nzchar(set_columns))) {
      stop("`set_columns` must be a non-empty character vector.", call. = FALSE)
    }
    missing_cols <- setdiff(set_columns, names(data))
    if (length(missing_cols)) {
      stop("`set_columns` missing column(s): ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    if (length(set_columns) < 2L || length(set_columns) > 3L) {
      stop("`plot_venn` supports two or three sets.", call. = FALSE)
    }
    sets <- lapply(set_columns, function(column) {
      idx <- !is.na(data[[column]]) & as.numeric(data[[column]]) == 1
      id[idx]
    })
    names(sets) <- set_columns
  } else if (is.list(data)) {
    if (!length(data)) stop("`data` list is empty.", call. = FALSE)
    if (length(data) < 2L || length(data) > 3L) {
      stop("`plot_venn` supports two or three sets.", call. = FALSE)
    }
    if (is.null(names(data)) || any(!nzchar(names(data)))) {
      stop("A list `data` must have non-empty set names.", call. = FALSE)
    }
    sets <- data
    set_columns <- names(sets)
  } else {
    stop("`data` must be a data frame or a named list.", call. = FALSE)
  }

  if (is.null(set_names)) {
    set_names <- set_columns
  }
  if (!is.character(set_names) || length(set_names) != length(sets) ||
      anyNA(set_names) || any(!nzchar(set_names))) {
    stop("`set_names` must match the number of sets.", call. = FALSE)
  }
  names(sets) <- set_names

  if (is.null(set_colors)) {
    palette <- .toolbox_group_palette()
    set_colors <- palette[seq_along(sets)]
  }
  if (!is.character(set_colors) || length(set_colors) != length(sets) ||
      anyNA(set_colors) || any(!nzchar(set_colors))) {
    stop("`set_colors` must match the number of sets.", call. = FALSE)
  }
  if (inherits(try(grDevices::col2rgb(set_colors), silent = TRUE), "try-error")) {
    stop("`set_colors` must contain valid R colors.", call. = FALSE)
  }

  plot <- ggvenn::ggvenn(
    sets,
    fill_color = set_colors,
    stroke_size = stroke_size,
    set_name_size = set_name_size,
    fill_alpha = fill_alpha,
    show_percentage = FALSE
  ) +
    ggplot2::labs(title = title)

  if (!is.null(style)) {
    if (!inherits(style, "bioinfo_plot_style")) {
      stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
    }
    plot <- plot + style$ggplot_theme
  } else {
    plot <- plot +
      ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

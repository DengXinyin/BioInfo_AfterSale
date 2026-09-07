# Read-QC, assembly and gene-length plots from the corresponding v2.87 scripts.

metage_base_content_plot <- function(data, position = "Position", title = "Base content",
                                     font_family = metage_font_family()) {
  if (!position %in% names(data)) stop("Position column is missing.")
  bases <- setdiff(names(data), position)
  if (!length(bases) || !all(vapply(data[bases], is.numeric, logical(1)))) {
    stop("Base columns must be numeric.")
  }
  dat <- data.frame(Position = rep(data[[position]], length(bases)),
                    Base = rep(bases, each = nrow(data)),
                    Percentage = unlist(data[bases], use.names = FALSE))
  ggplot2::ggplot(dat, ggplot2::aes(Position, Percentage, colour = Base)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(x = "Position along reads", y = "Percentage (%)", title = title) +
    .metage_theme(font_family)
}

metage_error_profile <- function(data, position = "Position", error = "ErrorRate",
                                  title = "Sequencing error rate",
                                  font_family = metage_font_family()) {
  if (!all(c(position, error) %in% names(data))) stop("Position/error columns are missing.")
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[position]], y = .data[[error]])) +
    ggplot2::geom_col(width = 0.8, fill = "#43CD80") +
    ggplot2::labs(x = "Position along reads", y = "Error rate", title = title) +
    .metage_theme(font_family)
}

metage_read_composition <- function(counts, relative = TRUE, colors = NULL,
                                     title = "Read quality composition",
                                     font_family = metage_font_family()) {
  x <- .metage_matrix(counts, "counts")
  if (relative) x <- sweep(x, 2, pmax(colSums(x), 1), "/")
  dat <- .metage_long(x, "Category", "Sample")
  dat$Category <- factor(dat$Category, levels = rownames(x))
  if (is.null(colors)) colors <- .metage_palette(nrow(x))
  ggplot2::ggplot(dat, ggplot2::aes(Sample, Value, fill = Category)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = if (relative) scales::percent else ggplot2::waiver(),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = if (relative) "Proportion" else "Read count", title = title) +
    .metage_theme(font_family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

metage_length_distribution <- function(lengths, breaks = NULL, group = NULL,
                                        title = "Sequence length distribution",
                                        font_family = metage_font_family()) {
  if (!is.numeric(lengths) || any(!is.finite(lengths))) stop("`lengths` must be finite numeric values.")
  if (is.null(breaks)) {
    dat <- data.frame(Length = lengths)
    return(ggplot2::ggplot(dat, ggplot2::aes(Length)) +
      ggplot2::geom_histogram(bins = 40, fill = "#377EB8", colour = "white") +
      ggplot2::labs(y = "Sequence count", title = title) + .metage_theme(font_family))
  }
  bins <- cut(lengths, breaks = breaks, include.lowest = TRUE, right = FALSE)
  dat <- data.frame(Bin = bins, Group = group %metage_or% "All")
  dat <- as.data.frame(table(dat), stringsAsFactors = FALSE)
  ggplot2::ggplot(dat, ggplot2::aes(Group, Freq, fill = Bin)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::labs(x = NULL, y = "Sequence count", title = title) +
    .metage_theme(font_family)
}

# Toolbox two-group survival palette (matches the 09_生存曲线 module).
.survival_group_palette <- function() {
  c("#E64B35", "#3182BD", "#00A087", "#F39B7F", "#7E6148")
}

#' Draw a Kaplan-Meier survival curve
#'
#' Fits a Kaplan-Meier curve with [survival::survfit()] and draws it with
#' [survminer::ggsurvplot()]. Supports confidence intervals, P-value/p-value
#' method labels, a risk table, and optional faceting by a covariate column.
#'
#' @param data Data frame containing the time, event, and group columns.
#' @param time_column Name of the time-to-event column.
#' @param event_column Name of the event indicator column (`1` event, `0` censor).
#' @param group_column Name of the grouping column used to split the strata.
#' @param palette Colour vector mapped to the group levels. `NULL` uses the
#'   survival palette.
#' @param conf_int Whether to draw the confidence interval bands.
#' @param pval, pval_method Whether to annotate the log-rank P value and the
#'   P-value method text.
#' @param risk_table Whether to draw the risk table below the curve.
#' @param facet_column Optional covariate column used with
#'   [survminer::ggsurvplot()]'s `facet.by`. `NULL` draws a single panel.
#' @param facet_ncol Number of facet columns when `facet_column` is set.
#' @param title Optional plot title.
#' @param x_label,y_label Axis labels.
#' @param legend_title Legend title.
#' @param font_family Font family passed to the ggsurvplot theme.
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return The [survminer::ggsurvplot()] list, whose `plot` element is the
#'   [ggplot2] object.
#' @export
plot_survival <- function(
    data,
    time_column = "time",
    event_column = "event",
    group_column = "group",
    palette = NULL,
    conf_int = TRUE,
    pval = TRUE,
    pval_method = TRUE,
    risk_table = TRUE,
    facet_column = NULL,
    facet_ncol = 3,
    title = NULL,
    x_label = "Time (Months)",
    y_label = "Survival Probability",
    legend_title = "Group",
    font_family = "Times New Roman",
    output_file = NULL,
    figure_width = 10,
    figure_height = 8,
    dpi = 300) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  for (column in c(time_column, event_column, group_column)) {
    if (!is.character(column) || length(column) != 1L || is.na(column) ||
        !nzchar(column) || !column %in% names(data)) {
      stop("Column `", column, "` is not a column in `data`.", call. = FALSE)
    }
  }
  if (!is.numeric(data[[time_column]])) {
    stop("`time_column` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(data[[event_column]])) {
    stop("`event_column` must be numeric.", call. = FALSE)
  }
  .assert_flag(conf_int, "conf_int")
  .assert_flag(pval, "pval")
  .assert_flag(pval_method, "pval_method")
  .assert_flag(risk_table, "risk_table")
  .assert_positive_number(facet_ncol, "facet_ncol")
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")
  facet_ncol <- as.integer(facet_ncol)
  if (facet_ncol < 1L) facet_ncol <- 1L

  group_levels <- unique(as.character(data[[group_column]]))
  if (is.null(palette)) {
    palette <- .survival_group_palette()
  }
  if (length(palette) < length(group_levels)) {
    palette <- rep(palette, length.out = length(group_levels))
  }

  formula <- stats::as.formula(
    paste("survival::Surv(", time_column, ",", event_column, ") ~", group_column)
  )
  fit <- eval(substitute(
    survival::survfit(formula, data = data), list(formula = formula)
  ))

  ggplot_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(text = ggplot2::element_text(family = font_family))

  facet_args <- list()
  if (!is.null(facet_column)) {
    if (!is.character(facet_column) || length(facet_column) != 1L ||
        !facet_column %in% names(data)) {
      stop("`facet_column` must be a column in `data`.", call. = FALSE)
    }
    facet_args$facet.by <- facet_column
    facet_args$ncol <- facet_ncol
  }

  ggsurv_args <- c(
    list(
      fit = fit,
      data = data,
      pval = pval,
      pval.method = pval_method,
      conf.int = conf_int,
      risk.table = risk_table,
      risk.table.col = "strata",
      palette = palette,
      xlab = x_label,
      ylab = y_label,
      title = title,
      legend.title = legend_title,
      ggtheme = ggplot_theme
    ),
    facet_args
  )

  surv_plot <- do.call(survminer::ggsurvplot, ggsurv_args)

  attr(surv_plot, "figure_width") <- figure_width
  attr(surv_plot, "figure_height") <- figure_height
  attr(surv_plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .survival_save_plot(surv_plot, output_file, figure_width, figure_height, dpi)
  }
  surv_plot
}

.survival_save_plot <- function(surv_plot, output_file, width, height, dpi) {
  if (!is.character(output_file) || length(output_file) != 1L ||
      is.na(output_file) || !nzchar(output_file)) {
    stop("`output_file` must be one non-empty filename.", call. = FALSE)
  }
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir, call. = FALSE)
  }
  extension <- tolower(tools::file_ext(output_file))
  if (!extension %in% c("pdf", "png")) {
    stop("`output_file` must end in .pdf or .png.", call. = FALSE)
  }
  if (extension == "pdf") {
    grDevices::cairo_pdf(
      output_file, width = width, height = height
    )
  } else {
    grDevices::png(
      output_file, width = width, height = height,
      units = "in", res = dpi, type = "cairo"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  print(surv_plot)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  invisible(NULL)
}

#' Draw a Cox regression forest plot
#'
#' Builds a forest plot of hazard ratios from a Cox proportional-hazards model
#' or its coefficient summary. When `model` is supplied, [survival::coxph()] is
#' assumed to have been run on `data` with `formula`; otherwise `hr_table` is
#' used directly and must contain `term`, `HR`, `lower`, and `upper` columns.
#'
#' @param hr_table Optional data frame of hazard-ratio results with columns
#'   `term`, `HR`, `lower`, `upper`, and optionally `p`. Provided when
#'   `model` is `NULL`.
#' @param model Optional [survival::coxph()] model object. When supplied, the
#'   table is computed from the model coefficients.
#' @param data Data frame used to fit `model` (required when `model` is given).
#' @param formula Cox model formula, e.g. `Surv(time, event) ~ group + age`.
#'   Required when `model` is given.
#' @param term_labels Optional named character vector mapping model term names
#'   to display labels.
#' @param title Optional plot title.
#' @param point_color Point colour for the hazard-ratio estimates.
#' @param x_label Axis label. Defaults to `"Hazard Ratio (95% CI)"`.
#' @param x_limits,x_breaks Axis limits and breaks (on the log10 scale).
#' @param font_family Font family used by the shared plot style. `"sans"` is the
#'   portable default.
#' @param style Optional object returned by [choose_plot_style()].
#' @param output_file Optional PDF or PNG output path.
#' @param figure_width,figure_height Output dimensions in inches.
#' @param dpi PNG resolution.
#'
#' @return A ggplot object.
#' @export
plot_forest <- function(
    hr_table = NULL,
    model = NULL,
    data = NULL,
    formula = NULL,
    term_labels = NULL,
    title = "Cox Regression Forest Plot",
    point_color = "#E64B35",
    x_label = "Hazard Ratio (95% CI)",
    x_limits = c(0.1, 30),
    x_breaks = c(0.1, 0.5, 1, 2, 5, 10, 20),
    font_family = "Times New Roman",
    style = NULL,
    output_file = NULL,
    figure_width = 10,
    figure_height = 6,
    dpi = 300) {

  if (!is.null(model)) {
    if (is.null(data) || is.null(formula)) {
      stop("`data` and `formula` are required when `model` is supplied.",
           call. = FALSE)
    }
    formula <- stats::as.formula(formula)
    table <- as.data.frame(stats::coef(summary(model)))
    table$term <- rownames(table)
    table$HR <- exp(table$coef)
    table$lower <- exp(table$coef - 1.96 * table[["se(coef)"]])
    table$upper <- exp(table$coef + 1.96 * table[["se(coef)"]])
    if ("Pr(>|z|)" %in% names(table)) table$p <- table[["Pr(>|z|)"]]
  } else {
    if (is.null(hr_table) || !is.data.frame(hr_table)) {
      stop("Provide either `hr_table` or `model`.", call. = FALSE)
    }
    for (column in c("term", "HR", "lower", "upper")) {
      if (!column %in% names(hr_table)) {
        stop("`hr_table` is missing column '", column, "'.", call. = FALSE)
      }
    }
    table <- as.data.frame(hr_table)
  }

  if (!is.character(point_color) || length(point_color) != 1L ||
      is.na(point_color) || !nzchar(point_color)) {
    stop("`point_color` must be one non-empty colour value.", call. = FALSE)
  }
  .assert_numeric_vector(x_limits, "x_limits", length_required = 2L, positive = TRUE)
  .assert_positive_number(figure_width, "figure_width")
  .assert_positive_number(figure_height, "figure_height")
  .assert_positive_number(dpi, "dpi")

  if (is.null(term_labels)) {
    table$label <- table$term
  } else {
    table$label <- unname(term_labels[table$term])
    table$label[is.na(table$label)] <- table$term[is.na(table$label)]
  }
  table$label <- factor(table$label, levels = rev(table$label))

  if (is.null(style)) {
    style <- choose_plot_style(
      font_family = font_family,
      title = list(size = 16),
      axis_title = list(size = 13),
      axis_text = list(size = 11),
      legend = list(show = FALSE)
    )
  }
  if (!inherits(style, "bioinfo_plot_style")) {
    stop("`style` must be created by `choose_plot_style()`.", call. = FALSE)
  }

  plot <- ggplot2::ggplot(
    table,
    ggplot2::aes(x = .data$HR, y = .data$label)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      orientation = "y", width = 0.2, color = "#333333"
    ) +
    ggplot2::geom_point(color = point_color, size = 3) +
    ggplot2::scale_x_log10(limits = x_limits, breaks = x_breaks) +
    ggplot2::labs(x = x_label, y = NULL, title = title) +
    style$ggplot_theme +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  attr(plot, "figure_width") <- figure_width
  attr(plot, "figure_height") <- figure_height
  attr(plot, "dpi") <- dpi

  if (!is.null(output_file)) {
    .pca_save_plot(plot, output_file, figure_width, figure_height, dpi)
  }
  plot
}

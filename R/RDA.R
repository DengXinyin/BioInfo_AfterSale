#' Plot a constrained RDA ordination
#' @param response Numeric response matrix/data frame with samples in rows.
#' @param environment Numeric explanatory-variable data frame with samples in rows.
#' @param group Optional sample grouping vector, named by sample ID when possible.
#' @param transform Response transformation: `"hellinger"` or `"none"`.
#' @param vif_cutoff Remove the highest VIF until all remaining VIFs are below this value.
#' @param max_variables Maximum number of explanatory variables retained.
#' @param permutations Number of permutations for the global test.
#' @param title Plot title.
#' @param style A style from [choose_plot_style()].
#' @param output_file Optional output filename.
#' @param width,height,dpi Optional output overrides.
#' @return A list containing plot, model, anova, selected_variables, site_scores, and environmental_vectors.
#' @export
plot_rda <- function(response, environment, group = NULL, transform = c("hellinger", "none"), vif_cutoff = 10, max_variables = 3, permutations = 999, title = "RDA", style = NULL, output_file = NULL, width = NULL, height = NULL, dpi = NULL) {
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Package `vegan` is required.", call. = FALSE)
  transform <- match.arg(transform); z <- .validate_style_output(style, output_file, width, height, dpi)
  Y <- as.data.frame(response, check.names = FALSE); E <- as.data.frame(environment, check.names = FALSE)
  if (is.null(rownames(Y)) || is.null(rownames(E))) stop("Both inputs must have sample row names.", call. = FALSE)
  common <- intersect(rownames(Y), rownames(E)); if (length(common) < 3L) stop("At least three shared samples are required.", call. = FALSE)
  Y <- Y[common, , drop = FALSE]; E <- E[common, , drop = FALSE]; Y[] <- lapply(Y, as.numeric); E[] <- lapply(E, as.numeric)
  if (any(!is.finite(as.matrix(Y))) || (transform == "hellinger" && any(as.matrix(Y) < 0))) stop(if (transform == "hellinger") "`response` must contain finite, non-negative values for Hellinger RDA." else "`response` must contain finite numeric values.", call. = FALSE)
  if (transform == "hellinger") Y <- vegan::decostand(Y, "hellinger")
  E <- E[, vapply(E, function(v) sd(v, na.rm = TRUE) > 0, logical(1)), drop = FALSE]; if (!ncol(E)) stop("No variable explanatory columns remain.", call. = FALSE)
  keep <- names(E); history <- data.frame(removed = character(), max_vif = numeric())
  repeat { fit0 <- vegan::rda(Y ~ ., data = E[, keep, drop = FALSE]); v <- tryCatch(vegan::vif.cca(fit0), error = function(e) rep(Inf, length(keep))); names(v) <- keep; if (length(keep) <= 1L || max(v, na.rm = TRUE) <= vif_cutoff) break; rmv <- names(which.max(v)); history <- rbind(history, data.frame(removed = rmv, max_vif = max(v, na.rm = TRUE))); keep <- setdiff(keep, rmv) }
  if (length(keep) > max_variables) keep <- keep[seq_len(max_variables)]
  fit <- vegan::rda(Y ~ ., data = E[, keep, drop = FALSE]); test <- vegan::anova.cca(fit, permutations = permutations)
  score_two <- function(display) { z <- tryCatch(as.data.frame(vegan::scores(fit, display = display, choices = 1:2, scaling = 2)), error = function(e) NULL); if (is.null(z)) z <- as.data.frame(vegan::scores(fit, display = display, choices = 1, scaling = 2)); if (ncol(z) == 1L) z$RDA2 <- 0; z[, seq_len(2), drop = FALSE] }
  site <- score_two("sites"); arrows <- score_two("bp"); site <- cbind(sample = rownames(site), site); names(site)[2:3] <- c("RDA1", "RDA2"); arrows <- cbind(variable = rownames(arrows), arrows); names(arrows)[2:3] <- c("RDA1", "RDA2")
  if (!is.null(group)) { if (!is.null(names(group))) group <- group[site$sample]; site$group <- factor(group) } else site$group <- factor("All samples")
  imp <- summary(fit)$concont$importance; pct <- as.numeric(imp[2, seq_len(min(2, ncol(imp)))]) * 100; pct <- c(pct, rep(NA_real_, 2 - length(pct))); pval <- as.numeric(test[1, "Pr(>F)"])
  p <- ggplot2::ggplot(site, ggplot2::aes(.data$RDA1, .data$RDA2, color = .data$group)) + ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") + ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey60") + ggplot2::geom_segment(data = arrows, ggplot2::aes(x = 0, y = 0, xend = .data$RDA1, yend = .data$RDA2), inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(.16, "cm")), color = "red") + ggplot2::geom_text(data = arrows, ggplot2::aes(x = .data$RDA1, y = .data$RDA2, label = .data$variable), inherit.aes = FALSE, vjust = -.5, color = "black") + ggplot2::geom_point(size = 3) + ggplot2::labs(title = title, subtitle = sprintf("Permutation P = %.4f", pval), x = sprintf("RDA1 (%.2f%%)", pct[1]), y = sprintf("RDA2 (%.2f%%)", pct[2]), color = "Group") + z$style$ggplot_theme
  .save_bioinfo_plot(p, z$style, output_file, z$width, z$height, z$dpi); list(plot = p, model = fit, anova = test, selected_variables = keep, vif_removal = history, site_scores = site, environmental_vectors = arrows)
}

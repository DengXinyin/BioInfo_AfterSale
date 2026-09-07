# Special plots consolidated from len_Dis.r and NCM.R in micro_v2:v5.42.

micro_length_distribution <- function(lengths, binwidth = 10,
                                      title = "Feature length distribution",
                                      font_family = "sans") {
  if (!is.numeric(lengths) || any(!is.finite(lengths))) stop("`lengths` must be finite numeric values.")
  dat <- data.frame(Length = lengths)
  ggplot2::ggplot(dat, ggplot2::aes(x = Length)) +
    ggplot2::geom_histogram(binwidth = binwidth, fill = "#377EB8", colour = "white") +
    ggplot2::labs(x = "Length (bp)", y = "Feature count", title = title) + .micro_theme(font_family)
}

micro_ncm_fit <- function(abundance) {
  x <- .micro_matrix(abundance)
  spp <- t(x)
  n_reads <- mean(rowSums(spp))
  p <- colMeans(spp) / n_reads
  freq <- colMeans(spp > 0)
  keep <- p > 0 & freq > 0 & p < 1
  p <- p[keep]; freq <- freq[keep]
  d <- 1 / n_reads
  objective <- function(log_m) {
    pred <- stats::pbeta(d, n_reads * exp(log_m) * p,
                         n_reads * exp(log_m) * (1 - p), lower.tail = FALSE)
    sum((freq - pred)^2)
  }
  opt <- stats::optimize(objective, interval = log(c(1e-8, 1e3)))
  m <- exp(opt$minimum)
  pred <- stats::pbeta(d, n_reads * m * p, n_reads * m * (1 - p), lower.tail = FALSE)
  se <- sqrt(pmax(pred * (1 - pred) / nrow(spp), 0))
  lower <- pmax(0, pred - 1.96 * se); upper <- pmin(1, pred + 1.96 * se)
  r2 <- 1 - sum((freq - pred)^2) / sum((freq - mean(freq))^2)
  axes <- data.frame(Feature = names(p), MeanRelativeAbundance = p,
                     Frequency = freq, Predicted = pred, Lower = lower, Upper = upper)
  structure(list(axes = axes, model = data.frame(R2 = r2, m = m, N = n_reads)),
            class = "micro_ncm")
}

micro_ncm_plot <- function(x, font_family = "sans") {
  fit <- if (inherits(x, "micro_ncm")) x else micro_ncm_fit(x)
  dat <- fit$axes[order(fit$axes$MeanRelativeAbundance), ]
  dat$Status <- ifelse(dat$Frequency < dat$Lower, "Below",
                       ifelse(dat$Frequency > dat$Upper, "Above", "Neutral"))
  ggplot2::ggplot(dat, ggplot2::aes(x = MeanRelativeAbundance, y = Frequency, colour = Status)) +
    ggplot2::geom_ribbon(
      data = dat,
      mapping = ggplot2::aes(x = MeanRelativeAbundance, ymin = Lower, ymax = Upper),
      inherit.aes = FALSE, fill = "#9ECAE1", alpha = 0.35
    ) +
    ggplot2::geom_line(ggplot2::aes(y = Predicted), colour = "#2166AC", linewidth = 1) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(Below = "#A52A2A", Neutral = "black", Above = "#29A6A6")) +
    ggplot2::labs(x = "Mean relative abundance (log10)", y = "Occurrence frequency",
                  title = "Neutral community model",
                  subtitle = sprintf("R2 = %.3f; Nm = %.1f", fit$model$R2, fit$model$m * fit$model$N)) +
    .micro_theme(font_family)
}

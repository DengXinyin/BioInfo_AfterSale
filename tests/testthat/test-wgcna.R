test_that("WGCNA expression preparation validates orientation and filtering", {
  skip_if_not_installed("WGCNA")
  set.seed(1)
  expression <- matrix(
    rexp(40 * 12), nrow = 40,
    dimnames = list(paste0("G", 1:40), paste0("S", 1:12))
  )
  expression[1, 1:5] <- 0
  prepared <- WGCNA_prepare_expression(expression, zero_fraction_max = 0.3)
  expect_s3_class(prepared, "WGCNA_prepared")
  expect_equal(dim(prepared$datExpr), c(12, 39))
  expect_equal(prepared$removed$zero_filter_genes, "G1")

  kept <- WGCNA_prepare_expression(expression, keep_all_genes = TRUE)
  expect_equal(ncol(kept$datExpr), 40)
  expect_equal(kept$parameters$transform, "none")
})

test_that("WGCNA image soft-power rule is explicit and reproducible", {
  skip_if_not_installed("WGCNA")
  set.seed(2)
  expression <- matrix(
    rnorm(80 * 15), nrow = 80,
    dimnames = list(paste0("G", 1:80), paste0("S", 1:15))
  )
  prepared <- WGCNA_prepare_expression(expression, zero_fraction_max = NULL)
  power <- WGCNA_select_power(
    prepared, powers = c(1, 3, 6, 9), strategy = "image_rule",
    network_type = "unsigned"
  )
  expect_s3_class(power, "WGCNA_power")
  expect_equal(power$power, 9)
  expect_match(power$reason, "reference image")
})

test_that("WGCNA modules, sample summaries, traits and kME stay aligned", {
  skip_if_not_installed("WGCNA")
  set.seed(3)
  n_samples <- 20
  latent_a <- rnorm(n_samples)
  latent_b <- rnorm(n_samples)
  expression <- rbind(
    t(replicate(30, latent_a + rnorm(n_samples, sd = 0.2))),
    t(replicate(30, latent_b + rnorm(n_samples, sd = 0.2)))
  )
  rownames(expression) <- paste0("G", seq_len(nrow(expression)))
  colnames(expression) <- paste0("S", seq_len(ncol(expression)))
  prepared <- WGCNA_prepare_expression(expression, zero_fraction_max = NULL)
  modules <- WGCNA_build_modules(
    prepared, power = 6, min_module_size = 10,
    max_block_size = 1000, threads = 1, verbose = 0
  )
  sample_info <- data.frame(
    Sample_ID = colnames(expression),
    time_group = rep(c("A", "B"), each = 10),
    stringsAsFactors = FALSE
  )
  sample_result <- WGCNA_module_sample(modules, sample_info, "time_group")
  membership <- WGCNA_module_membership(modules, sample_result$MEs_oriented)
  traits <- data.frame(Sample_ID = colnames(expression), response = latent_a)
  trait_result <- WGCNA_module_trait(
    modules, traits, MEs = sample_result$MEs_oriented
  )

  expect_equal(nrow(modules$module_gene), 60)
  expect_equal(nrow(sample_result$sample_values),
               n_samples * ncol(sample_result$MEs_oriented))
  expect_equal(nrow(sample_result$time_summary),
               2 * ncol(sample_result$MEs_oriented))
  expect_equal(nrow(membership$assigned), 60)
  expect_equal(ncol(trait_result$correlation), 1)
})

test_that("WGCNA no-trait workflow writes reusable core tables", {
  skip_if_not_installed("WGCNA")
  set.seed(4)
  expression <- matrix(
    rexp(60 * 12), nrow = 60,
    dimnames = list(paste0("G", 1:60), paste0("S", 1:12))
  )
  output <- tempfile("wgcna-output-")
  result <- WGCNA_run(
    expression,
    prepare_args = list(zero_fraction_max = NULL),
    power_args = list(strategy = "manual", manual_power = 6, powers = c(1, 6)),
    module_args = list(min_module_size = 10, threads = 1),
    output_dir = output,
    verbose = FALSE
  )
  expect_s3_class(result, "WGCNA_workflow")
  expect_true(all(file.exists(file.path(
    output,
    c("expression_qc.csv", "soft_threshold_fit.csv", "modulegene.csv",
      "module_sample_eigengene_long.csv", "assigned_module_kME.csv",
      "RUN_MANIFEST.csv")
  ))))
})

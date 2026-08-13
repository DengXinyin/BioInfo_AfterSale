test_that("unified interface requires exactly one input", {
  expect_error(enrichment_plot(), "exactly one")
  expect_error(
    enrichment_plot(gene = "1", result = data.frame()),
    "exactly one"
  )
})

test_that("run_enrichment validates genes before loading annotation", {
  expect_error(run_enrichment(character(), analysis = "GO"), "no usable")
})

test_that("plot_enrichment validates the result type", {
  expect_error(plot_enrichment(data.frame()), "enrichResult")
})

test_that("unified interface requires exactly one input", {
  expect_error(GO_KEGG(), "exactly one")
  expect_error(
    GO_KEGG(gene = "1", result = data.frame()),
    "exactly one"
  )
})

test_that("GO_KEGG_analyse validates genes before loading annotation", {
  expect_error(GO_KEGG_analyse(character(), analysis = "GO"), "no usable")
})

test_that("GO_KEGG_plot validates the result type", {
  expect_error(GO_KEGG_plot(data.frame()), "enrichResult")
})

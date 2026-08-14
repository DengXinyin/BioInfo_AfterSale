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

test_that("choose_plot_style creates and validates a shared style", {
  style <- choose_plot_style(
    font_family = "Arial",
    title = list(size = 22, bold = TRUE),
    axis_text = list(size = 15),
    legend = list(position = "bottom")
  )
  expect_s3_class(style, "bioinfo_plot_style")
  expect_s3_class(style$ggplot_theme, "theme")
  expect_equal(style$global$font_family, "Arial")
  expect_equal(style$text$axis_text$size, 15)
  expect_equal(style$legend$position, "bottom")
  expect_error(choose_plot_style(theme = "invalid"))
  expect_error(choose_plot_style(title = list(unknown = 1)), "Unknown")
})

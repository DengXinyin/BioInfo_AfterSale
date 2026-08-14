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

test_that("GO_KEGG_plot validates a data-frame result", {
  expect_error(GO_KEGG_plot(data.frame()), "missing column")
  expect_error(GO_KEGG_plot(1), "enrichResult")
})

test_that("GO_KEGG_plot supports custom table mappings", {
  table <- data.frame(
    Description = c("Pathway A", "Pathway B"),
    pvalue = c(0.001, 0.02),
    RichFactor = c(2.5, 4.2),
    Count = c(5, 9)
  )
  style <- choose_plot_style(font_family = "sans")
  plot <- GO_KEGG_plot(
    table,
    filter_by = NULL,
    x = "pvalue",
    x_transform = "neg_log10",
    color = "RichFactor",
    style = style
  )
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 2)
  expect_equal(sort(plot$data$.XValue), sort(-log10(table$pvalue)))
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

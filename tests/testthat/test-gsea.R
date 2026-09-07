test_that("GSEA running score is reproducible and ends at zero", {
  ranked <- c(G1 = 3, G2 = 2, G3 = 1, G4 = -1, G5 = -2, G6 = -3)
  curve <- calculate_gsea_curve(ranked, c("G1", "G3"), exponent = 1)

  expect_s3_class(curve, "bioinfo_gsea_curve")
  expect_equal(curve$hits, c(1L, 3L))
  expect_equal(tail(curve$running_score, 1), 0)
  expect_equal(curve$peak_position, which.max(abs(curve$running_score)))
  expect_true(all(curve$leading_edge %in% c("G1", "G3")))
})

test_that("GSEA running score validates gene identifiers and overlap", {
  ranked <- setNames(c(2, 1, -1), c("A", "B", "C"))
  expect_error(calculate_gsea_curve(unname(ranked), "A"), "named numeric")
  expect_error(calculate_gsea_curve(setNames(ranked, c("A", "A", "C")), "A"), "unique")
  expect_error(calculate_gsea_curve(ranked, "missing"), "no member")
  expect_error(calculate_gsea_curve(ranked, names(ranked)), "every ranked gene")
  expect_error(calculate_gsea_curve(ranked, "A", exponent = -1), "non-negative")
})

test_that("GSEA pathway wrapper adds metric band, peak, and curve metadata", {
  set.seed(10)
  ranked <- sort(setNames(rnorm(200), paste0("G", 1:200)), decreasing = TRUE)
  genes <- names(ranked)[c(1:15, 80:90)]
  plot <- plot_gsea_pathway(
    ranked, genes,
    statistics = c(NES = 1.8, `P value` = 0.002, `Adjusted P` = 0.01),
    style = choose_plot_style(font_family = "sans")
  )

  expect_s3_class(plot, "ggplot")
  expect_s3_class(attr(plot, "gsea_curve"), "bioinfo_gsea_curve")
  expect_equal(attr(plot, "peak_position"), attr(plot, "gsea_curve")$peak_position)
  expect_gte(length(plot$layers), 8)
})

test_that("low-level GSEA plot validates logical hits and named statistics", {
  score <- seq(0, 1, length.out = 10)
  metric <- seq(2, -2, length.out = 10)
  style <- choose_plot_style(font_family = "sans")
  expect_error(plot_gsea(score, rep(TRUE, 3), metric, style = style), "Logical")
  expect_error(plot_gsea(score, 2:3, metric, statistics = 1, style = style), "named")
})

test_that("gseaResult wrappers extract pathways and batch-export figures", {
  if (!methods::isClass("gseaResult")) {
    methods::setClass(
      "gseaResult",
      slots = c(result = "data.frame", geneList = "numeric", geneSets = "list")
    )
  }
  ranked <- sort(setNames(seq(3, -3, length.out = 100), paste0("G", 1:100)),
                 decreasing = TRUE)
  result_table <- data.frame(
    ID = "path:001", Description = "Synthetic pathway",
    NES = 1.75, pvalue = 0.002, p.adjust = 0.01,
    stringsAsFactors = FALSE
  )
  object <- methods::new("gseaResult")
  methods::slot(object, "result") <- result_table
  methods::slot(object, "geneList") <- ranked
  methods::slot(object, "geneSets") <- list(`path:001` = names(ranked)[1:20])

  plot <- plot_gsea_result(
    object, "path:001", style = choose_plot_style(font_family = "sans")
  )
  expect_s3_class(plot, "ggplot")
  expect_equal(attr(plot, "gsea_curve")$gene_set, names(ranked)[1:20])

  output <- tempfile("gsea-batch-")
  files <- plot_gsea_batch(
    object, "path:001", output, style = choose_plot_style(font_family = "sans")
  )
  expect_length(files, 1)
  expect_true(file.exists(files))
  expect_gt(file.info(files)$size, 0)
})

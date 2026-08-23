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

test_that("GO_KEGG_plot keeps only positive P values below the cutoff", {
  table <- data.frame(
    Description = c("Underflow", "Significant", "At cutoff", "Not significant"),
    pvalue = c(0, 0.01, 0.05, 0.10),
    RichFactor = c(0.1, 0.2, 0.3, 0.4),
    Count = c(1, 2, 3, 4)
  )
  expect_warning(
    plot <- GO_KEGG_plot(
      table, filter_by = "pvalue", cutoff = 0.05,
      x = "RichFactor", color = "pvalue", font_family = "sans"
    ),
    "pvalue = 0"
  )
  expect_equal(nrow(plot$data), 1)
  expect_identical(as.character(plot$data$.Label), "Significant")
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
  expect_equal(plot$scales$get_scales("size")$breaks, c(5, 9))
})

test_that("GO_KEGG_plot exposes after-sales axis, point, and guide controls", {
  table <- data.frame(
    Description = c("Pathway A", "Pathway B"),
    pvalue = c(0.001, 0.02),
    RichFactor = c(2.5, 4.2),
    Count = c(5, 9)
  )
  plot <- GO_KEGG_plot(
    table,
    filter_by = NULL,
    x = "pvalue",
    x_transform = "neg_log10",
    color = "RichFactor",
    x_label = expression(-log[10](Pvalue)),
    y_label = "Selected pathway",
    size_breaks = c(5, 7, 9),
    size_range = c(2, 8),
    point_alpha = 0.75,
    x_limits = c(1, 4),
    x_breaks = 1:4,
    x_expand = c(0, 0),
    legend_order = c("color", "size"),
    style = choose_plot_style(font_family = "sans")
  )
  size_scale <- plot$scales$get_scales("size")
  x_scale <- plot$scales$get_scales("x")
  expect_equal(plot$labels$y, "Selected pathway")
  expect_equal(size_scale$breaks, c(5, 7, 9))
  expect_equal(size_scale$palette(c(0, 1)), c(2, 8))
  expect_equal(plot$layers[[1]]$aes_params$alpha, 0.75)
  expect_equal(x_scale$limits, c(1, 4))
  expect_equal(x_scale$breaks, 1:4)
  expect_equal(x_scale$expand, c(0, 0))
  expect_error(GO_KEGG_plot(table, size_range = c(8, 2)), "minimum")
  expect_error(GO_KEGG_plot(table, legend_order = c("size", "size")), "unique")
})

test_that("GO_KEGG_plot defaults to at most three Gene count breaks", {
  table <- data.frame(
    Description = paste("Pathway", LETTERS[1:6]),
    pvalue = c(0.001, 0.002, 0.004, 0.008, 0.016, 0.032),
    RichFactor = c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45),
    Count = c(5, 7, 9, 11, 13, 15)
  )
  plot <- GO_KEGG_plot(
    table,
    filter_by = NULL,
    x = "pvalue",
    x_transform = "neg_log10",
    color = "RichFactor",
    font_family = "sans"
  )
  expect_equal(plot$scales$get_scales("size")$breaks, c(5, 9, 15))
})

test_that("GO_KEGG_plot highlights selected term labels", {
  table <- data.frame(
    Description = c("DNA repair", "Cell cycle", "RNA processing"),
    pvalue = c(0.001, 0.01, 0.02),
    RichFactor = c(0.42, 0.35, 0.28),
    Count = c(12, 9, 6)
  )
  plot <- GO_KEGG_plot(
    table,
    filter_by = NULL,
    x = "RichFactor",
    color = "pvalue",
    highlight_terms = c("DNA repair", "RNA processing"),
    highlight_color = "#7B2CBF",
    highlight_bold = TRUE,
    font_family = "sans"
  )
  labels <- as.character(plot$data$.Label)
  expect_true(any(grepl("DNA repair", labels, fixed = TRUE)))
  expect_true(any(grepl("RNA processing", labels, fixed = TRUE)))
  expect_equal(sum(grepl("font-weight:bold", labels, fixed = TRUE)), 2)
  expect_equal(sum(grepl("color:#7B2CBF", labels, fixed = TRUE)), 2)
  expect_s3_class(plot$theme$axis.text.y, "element_markdown")
  expect_warning(
    GO_KEGG_plot(
      table,
      filter_by = NULL,
      x = "RichFactor",
      color = "pvalue",
      highlight_terms = "Missing pathway",
      font_family = "sans"
    ),
    "not displayed"
  )
  expect_error(
    GO_KEGG_plot(
      table,
      filter_by = NULL,
      x = "RichFactor",
      color = "pvalue",
      highlight_terms = "DNA repair",
      highlight_color = "not-a-color"
    ),
    "valid R color"
  )
})

test_that("GO_KEGG_plot applies count breaks and highlights to enrichResult", {
  result_table <- data.frame(
    ID = paste0("GO:", 1:4),
    Description = c("DNA repair", "Cell cycle", "RNA processing", "Oxidative stress"),
    GeneRatio = c("5/20", "4/20", "3/20", "2/20"),
    BgRatio = rep("10/100", 4),
    RichFactor = c(0.50, 0.40, 0.30, 0.20),
    FoldEnrichment = c(2.5, 2.0, 1.5, 1.0),
    zScore = rep(1, 4),
    pvalue = c(0.001, 0.005, 0.010, 0.020),
    p.adjust = c(0.004, 0.010, 0.020, 0.030),
    qvalue = c(0.004, 0.010, 0.020, 0.030),
    geneID = c("1/2/3/4/5", "1/2/3/4", "1/2/3", "1/2"),
    Count = c(5, 4, 3, 2),
    stringsAsFactors = FALSE
  )
  enrichment <- methods::new(
    "enrichResult",
    result = result_table,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    organism = "human",
    ontology = "BP",
    gene = as.character(1:5),
    keytype = "ENTREZID",
    universe = as.character(1:100),
    geneSets = list(),
    readable = FALSE
  )
  plot <- GO_KEGG_plot(
    enrichment,
    filter_by = NULL,
    show_category = 4,
    font_family = "sans",
    highlight_terms = "DNA repair",
    highlight_color = "#7B2CBF"
  )
  expect_equal(plot$scales$get_scales("size")$breaks, c(2, 3, 5))
  expect_s3_class(plot$theme$axis.text.y, "element_markdown")
  y_scale <- plot$scales$get_scales("y")
  formatted <- y_scale$labels(c("DNA repair", "Cell cycle"))
  expect_match(formatted[[1]], "font-weight:bold")
  expect_match(formatted[[1]], "color:#7B2CBF")
  expect_identical(formatted[[2]], "Cell cycle")
})

test_that("GO_KEGG_plot saves selectable inch dimensions", {
  table <- data.frame(
    Description = c("Pathway A", "Pathway B"),
    pvalue = c(0.001, 0.02),
    RichFactor = c(2.5, 4.2),
    Count = c(5, 9)
  )
  output <- tempfile(fileext = ".pdf")
  style <- choose_plot_style(font_family = "sans")
  plot <- GO_KEGG_plot(
    table,
    filter_by = NULL,
    x = "pvalue",
    x_transform = "neg_log10",
    color = "RichFactor",
    style = style,
    output_file = output,
    figure_width = 11,
    figure_height = 6,
    dpi = 200
  )
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_equal(attr(plot, "figure_width"), 11)
  expect_equal(attr(plot, "figure_height"), 6)
  expect_equal(attr(plot, "dpi"), 200)
})

test_that("choose_plot_style creates and validates a shared style", {
  style <- choose_plot_style(
    font_family = "Arial",
    title = list(size = 22, bold = TRUE),
    axis_text = list(size = 15),
    legend = list(position = "bottom", box = "horizontal"),
    panel = list(
      border = TRUE,
      major_grid = TRUE,
      minor_grid = TRUE,
      border_width = 1,
      major_grid_color = "grey70"
    )
  )
  expect_s3_class(style, "bioinfo_plot_style")
  expect_s3_class(style$ggplot_theme, "theme")
  expect_equal(style$global$font_family, "Arial")
  expect_equal(style$text$axis_text$size, 15)
  expect_equal(style$legend$position, "bottom")
  expect_equal(style$legend$box, "horizontal")
  expect_true(style$panel$border)
  expect_true(style$panel$major_grid)
  expect_true(style$panel$minor_grid)
  expect_s3_class(style$ggplot_theme$panel.border, "element_rect")
  expect_s3_class(style$ggplot_theme$panel.grid.major, "element_line")
  expect_error(choose_plot_style(theme = "invalid"))
  expect_error(choose_plot_style(title = list(unknown = 1)), "Unknown")
  expect_error(choose_plot_style(panel = list(border = "yes")), "TRUE or FALSE")
})

test_that("Heatmap_plot exposes clustering, labels, groups, legends, and titles", {
  matrix <- matrix(
    c(1, 2, 3, 4, 5, 6, 6, 5, 4, 3, 2, 1),
    nrow = 3,
    dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("S", 1:4))
  )
  group <- c(S1 = "Control", S2 = "Control", S3 = "Treatment", S4 = "Treatment")
  heatmap <- Heatmap_plot(
    matrix,
    group = group,
    group_colors = c(Control = "grey50", Treatment = "brown"),
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_font_family = "sans",
    column_names_font_family = "sans",
    row_names_italic = TRUE,
    column_names_rot = 45,
    title = "Selected genes",
    title_position = "top",
    title_font_family = "sans",
    title_font_size = 16,
    title_font_face = "bold",
    legend_side = "right",
    draw_plot = FALSE
  )
  expect_s4_class(heatmap, "Heatmap")
  expect_true(heatmap@row_names_param$show)
  expect_true(heatmap@column_names_param$show)
  expect_false(heatmap@column_dend_param$cluster)
  expect_equal(heatmap@column_names_param$rot, 45)
  expect_equal(heatmap@column_title, "Selected genes")
  expect_s4_class(heatmap@top_annotation, "HeatmapAnnotation")
  expect_equal(attr(heatmap, "aftersale_draw_args")$heatmap_legend_side, "right")
  expect_equal(attr(heatmap, "aftersale_draw_args")$annotation_legend_side, "right")
  expect_true(attr(heatmap, "aftersale_draw_args")$merge_legends)
})

test_that("Heatmap_plot saves PDF output", {
  matrix <- matrix(rnorm(24), nrow = 6)
  rownames(matrix) <- paste0("Gene", 1:6)
  colnames(matrix) <- paste0("Sample", 1:4)
  output <- tempfile(fileext = ".pdf")
  Heatmap_plot(
    matrix,
    group = c("A", "A", "B", "B"),
    row_names_font_family = "sans",
    column_names_font_family = "sans",
    title_font_family = "sans",
    output_file = output,
    figure_width = 6,
    figure_height = 5,
    draw_plot = FALSE
  )
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  expect_error(
    Heatmap_plot(matrix, group = c("A", "B"), draw_plot = FALSE),
    "one value per matrix column"
  )
  expect_error(
    Heatmap_plot(matrix, zscore_breaks = c(0, -1, 1), draw_plot = FALSE),
    "strictly increasing"
  )
})

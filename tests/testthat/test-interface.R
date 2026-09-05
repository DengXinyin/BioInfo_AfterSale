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

test_that("volcano_plot classifies genes and adds three threshold lines", {
  table <- data.frame(
    log2FoldChange = c(-2, -1, 0, 1, 2, NA),
    padj = c(0.01, 0.01, 0.01, 0.05, 0.001, 0)
  )
  expect_warning(
    plot <- volcano_plot(table, font_family = "sans"),
    "non-positive padj"
  )
  expect_equal(table(plot$data$.Group), c(Down = 1, `Not significant` = 1, Up = 2))
  expect_equal(length(plot$layers), 3)
  expect_equal(plot$scales$get_scales("colour")$breaks,
               c("Up", "Down", "Not significant"))
  expect_equal(attr(plot, "log2fc_cutoff"), 1)
  expect_equal(attr(plot, "padj_cutoff"), 0.05)
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

test_that("plot_pca computes scores and stores group/sample columns", {
  set.seed(42)
  ng <- 20; ns <- 60
  expr <- matrix(rnorm(ng * ns, mean = 5, sd = 1.5), nrow = ng, ncol = ns)
  group <- rep(c("Control", "Treatment_A", "Treatment_B"), each = 20)
  expr[1:4, group == "Treatment_A"] <- expr[1:4, group == "Treatment_A"] + 2
  expr[1:4, group == "Treatment_B"] <- expr[1:4, group == "Treatment_B"] - 2
  colnames(expr) <- paste0("S", seq_len(ns))
  rownames(expr) <- paste0("g", seq_len(ng))
  df <- as.data.frame(t(expr))
  df$Sample <- colnames(expr)
  df$Group <- group

  plot <- plot_pca(df, sample_column = "Sample", group_column = "Group",
                   title = "PCA pilot")
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), ns)
  expect_true(all(c(".PC1", ".PC2", ".Group", ".Sample") %in% names(plot$data)))
  expect_equal(levels(plot$data$.Group),
               c("Control", "Treatment_A", "Treatment_B"))
  expect_equal(attr(plot, "x_pc"), 1)
  expect_equal(attr(plot, "y_pc"), 2)
  expect_match(plot$labels$x, "PC1 \\([0-9.]+%\\)")
  expect_match(plot$labels$y, "PC2 \\([0-9.]+%\\)")
  colour_scale <- plot$scales$get_scales("colour")
  expect_equal(colour_scale$name, "Group")
  expect_equal(colour_scale$breaks, c("Control", "Treatment_A", "Treatment_B"))
  # default toolbox palette is applied per group (visible in built data)
  built <- ggplot2::ggplot_build(plot)
  point_data <- built$data[[1]]
  expect_setequal(unique(point_data$colour),
                  c("#4DBBD5", "#E64B35", "#00A087"))
})

test_that("plot_pca supports custom PC pair and sample labels", {
  set.seed(1)
  ng <- 15; ns <- 30
  expr <- matrix(rnorm(ng * ns), nrow = ng, ncol = ns)
  group <- rep(c("A", "B"), each = 15)
  colnames(expr) <- paste0("s", seq_len(ns))
  rownames(expr) <- paste0("g", seq_len(ng))
  df <- as.data.frame(t(expr))
  df$Sample <- colnames(expr)
  df$Group <- group

  plot <- plot_pca(df, x_pc = 3, y_pc = 4, show_labels = TRUE,
                   group_colors = c(A = "#E64B35", B = "#4DBBD5"))
  expect_equal(attr(plot, "x_pc"), 3)
  expect_equal(attr(plot, "y_pc"), 4)
  expect_match(plot$labels$x, "PC3 \\([0-9.]+%\\)")
  expect_match(plot$labels$y, "PC4 \\([0-9.]+%\\)")
  expect_true(any(vapply(plot$layers, function(l) inherits(l$geom, "GeomText"),
                         logical(1))))
  built <- ggplot2::ggplot_build(plot)
  expect_setequal(unique(built$data[[1]]$colour), c("#E64B35", "#4DBBD5"))
})

test_that("plot_pca validates inputs and reports missing columns", {
  set.seed(7)
  df <- data.frame(Sample = paste0("S", 1:6),
                   Group = rep(c("A", "B"), each = 3),
                   x = rnorm(6), y = rnorm(6), z = rnorm(6))
  expect_error(plot_pca(df, group_column = "Nope"), "missing column")
  expect_error(plot_pca(df, expr_columns = c("x", "nope")), "missing column")
  expect_error(plot_pca(df, expr_columns = "x"), "At least two")
  expect_error(plot_pca(df, point_size = -1), "positive")
  expect_error(plot_pca(df, x_pc = 1.5), "whole number")
})
test_that("plot_venn draws a three-set Venn from indicator columns", {
  df <- data.frame(
    Gene = paste0("g", 1:12),
    DEG_SetA = c(1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
    DEG_SetB = c(0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0),
    DEG_SetC = c(0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0)
  )
  plot <- plot_venn(df, id_column = "Gene", title = "Venn")
  expect_s3_class(plot, "ggplot")
  expect_equal(attr(plot, "figure_width"), 8)
  plot2 <- plot_venn(df, id_column = "Gene",
                     set_colors = c("#E64B35", "#4DBBD5", "#00A087"))
  expect_s3_class(plot2, "ggplot")
})

test_that("plot_venn accepts a list input and renames sets", {
  sets <- list(A = letters[1:5], B = letters[3:7], C = letters[6:9])
  plot <- plot_venn(sets, set_names = c("Set A", "Set B", "Set C"))
  expect_s3_class(plot, "ggplot")
  expect_error(plot_venn(list(A = letters[1:3])), "two or three")
})

test_that("plot_venn validates inputs", {
  df <- data.frame(id = paste0("g", 1:6), a = c(1, 1, 0, 0, 0, 0), b = c(0, 0, 1, 1, 0, 0))
  expect_error(plot_venn(df, id_column = "nope"), "not a column")
  expect_error(plot_venn(df, id_column = "id", set_columns = c("x", "y")),
               "missing column")
  expect_error(plot_venn(df, id_column = "id", set_columns = c("a", "b"),
                         set_colors = c("red", "blue", "green")),
               "must match the number of sets")
  expect_error(plot_venn(df, fill_alpha = 2), "between 0 and 1")
  expect_error(plot_venn(1), "data frame or a named list")
})
test_that("plot_violin_box draws a violin+box+jitter plot", {
  df <- data.frame(
    Species = rep(c("setosa", "versicolor", "virginica"), each = 10),
    Feature = rep("Sepal.Length", 30),
    Value = c(rnorm(10, 5), rnorm(10, 6), rnorm(10, 6.5))
  )
  plot <- plot_violin_box(df, value_column = "Value", group_column = "Species",
                          title = "Violin + Box + Jitter")
  expect_s3_class(plot, "ggplot")
  labs <- ggplot2::get_labs(plot)
  expect_equal(labs$x, "Species")
  expect_equal(labs$y, "Value")
  expect_equal(labs$title, "Violin + Box + Jitter")
  # renders to a file without error (accepts a real data path)
  tmp <- tempfile(fileext = ".pdf")
  plot_violin_box(df, value_column = "Value", group_column = "Species",
                  output_file = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_violin_box supports faceting and plot type", {
  df <- data.frame(
    Species = rep(c("setosa", "virginica"), each = 20),
    Feature = rep(c("Sepal.Length", "Petal.Length"), each = 20),
    Value = rnorm(40)
  )
  plot <- plot_violin_box(df, value_column = "Value", group_column = "Species",
                          facet_column = "Feature", plot_type = "violin",
                          show_jitter = FALSE)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$facet, "FacetWrap")
  b <- ggplot2::ggplot_build(plot)
  expect_true(length(b$data) >= 1)
})

test_that("plot_violin_box validates inputs", {
  df <- data.frame(Species = letters[1:6], x = letters[1:6], Value = rnorm(6))
  expect_error(plot_violin_box(df, value_column = "nope"), "not a column")
  expect_error(plot_violin_box(df, value_column = "x"), "must be numeric")
  expect_error(plot_violin_box(df, group_column = "nope"), "not a column")
  expect_error(plot_violin_box(df, jitter_alpha = -1), "between 0 and 1")
})
test_that("plot_scatter draws an overall scatter with correlation annotation", {
  df <- data.frame(VarX = rnorm(30), VarY = rnorm(30) + 0.5)
  plot <- plot_scatter(df, x_column = "VarX", y_column = "VarY",
                       title = "Overall")
  expect_s3_class(plot, "ggplot")
  labs <- ggplot2::get_labs(plot)
  expect_equal(labs$x, "VarX")
  expect_equal(labs$y, "VarY")
  expect_equal(labs$title, "Overall")
  tmp <- tempfile(fileext = ".pdf")
  plot_scatter(df, x_column = "VarX", y_column = "VarY", output_file = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_scatter supports grouping and custom correlations", {
  df <- data.frame(Group = rep(c("A", "B"), each = 20),
                   VarX = rnorm(40), VarY = rnorm(40))
  plot <- plot_scatter(df, x_column = "VarX", y_column = "VarY",
                       group_column = "Group", cor_method = "spearman",
                       smooth_method = "loess")
  expect_s3_class(plot, "ggplot")
  expect_equal(ggplot2::get_labs(plot)$title, NULL)
  b <- ggplot2::ggplot_build(plot)
  expect_true(length(b$data) >= 2)
})

test_that("plot_scatter supports faceting", {
  df <- data.frame(Group = rep(c("A", "B", "C"), each = 15),
                   VarX = rnorm(45), VarY = rnorm(45))
  plot <- plot_scatter(df, x_column = "VarX", y_column = "VarY",
                       facet_column = "Group")
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$facet, "FacetWrap")
})

test_that("plot_scatter validates inputs", {
  df <- data.frame(VarX = rnorm(6), VarY = rnorm(6))
  expect_error(plot_scatter(df, x_column = "nope"), "not a column")
  expect_error(plot_scatter(df, point_alpha = -1), "between 0 and 1")
  expect_error(plot_scatter(df, point_size = 0), "positive")
  expect_scatter_missing <- tryCatch({ plot_scatter(df, group_column = "nope") },
                                     error = function(e) e)
  expect_true(inherits(expect_scatter_missing, "error"))
})
test_that("plot_distribution draws histogram, density, and both", {
  df <- data.frame(Expression = rnorm(100, 5, 1), Group = rep(c("Control", "Treatment"), each = 50))
  ph <- plot_distribution(df, plot_type = "histogram", title = "Hist")
  expect_s3_class(ph, "ggplot")
  expect_equal(ggplot2::get_labs(ph)$title, "Hist")
  pd <- plot_distribution(df, plot_type = "density")
  expect_s3_class(pd, "ggplot")
  pb <- plot_distribution(df, plot_type = "both")
  expect_s3_class(pb, "ggplot")
  tmp <- tempfile(fileext = ".pdf")
  plot_distribution(df, output_file = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_distribution supports grouping by a second column", {
  df <- data.frame(Expression = c(rnorm(50, 4), rnorm(50, 6)),
                   Group = rep(c("Control", "Treatment"), each = 50))
  plot <- plot_distribution(df, group_column = "Group")
  expect_s3_class(plot, "ggplot")
  b <- ggplot2::ggplot_build(plot)
  expect_true(length(b$data) >= 1)
})

test_that("plot_distribution validates inputs", {
  df <- data.frame(Expression = rnorm(20), Group = rep(c("A", "B"), 10))
  expect_error(plot_distribution(df, value_column = "nope"), "not a column")
  expect_error(plot_distribution(df, value_column = "Group"), "must be numeric")
  expect_error(plot_distribution(df, bins = 0), "positive")
  expect_error(plot_distribution(df, group_column = "nope"), "must be a column")
  expect_error(plot_distribution(df, fill_color = ""), "colour")
})
test_that("plot_survival fits a Kaplan-Meier curve", {
  df <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8),
    event = c(1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1),
    group = rep(c("Treat", "Placebo"), each = 8)
  )
  res <- plot_survival(df, time_column = "time", event_column = "event",
                       group_column = "group", title = "KM")
  expect_s3_class(res, "ggsurvplot")
  expect_true(inherits(res[["plot"]], "ggplot"))
  expect_equal(attr(res, "figure_width"), 10)
})

test_that("plot_survival supports faceting", {
  df <- data.frame(
    time = rep(c(1, 2, 3, 4, 5, 6, 7, 8), 4),
    event = rep(c(1, 0, 1, 1, 0, 1, 1, 0), 4),
    group = rep(c("A", "B"), each = 16),
    stage = rep(c("I", "II", "I", "II"), each = 8)
  )
  res <- plot_survival(df, group_column = "group", facet_column = "stage")
  # faceting makes ggsurvplot return a ggplt object
  expect_true(is(res, "ggplot") || is(res, "ggsurvplot"))
})

test_that("plot_forest builds a Cox forest plot from a model", {
  df <- data.frame(
    time = c(4, 6, 5, 8, 3, 7, 9, 2),
    event = c(1, 1, 0, 1, 1, 0, 1, 1),
    group = c("A", "B", "A", "B", "A", "B", "A", "B"),
    age = c(60, 55, 70, 65, 50, 68, 72, 58)
  )
  model <- survival::coxph(survival::Surv(time, event) ~ group + age, data = df)
  plot <- plot_forest(model = model, data = df,
                      formula = survival::Surv(time, event) ~ group + age)
  expect_s3_class(plot, "ggplot")
  expect_equal(ggplot2::get_labs(plot)$title, "Cox Regression Forest Plot")
})

test_that("plot_forest works from a precomputed hr_table", {
  hr_table <- data.frame(
    term = c("groupB", "age"),
    HR = c(1.5, 1.1), lower = c(0.9, 0.9), upper = c(2.5, 1.4)
  )
  plot <- plot_forest(hr_table = hr_table)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_survival and plot_forest validate inputs", {
  expect_error(plot_survival(data.frame(), group_column = "g"),
               "not a column")
  expect_error(plot_forest(), "Provide either")
  expect_error(plot_forest(hr_table = data.frame(term = "x")), "missing column")
})
test_that("plot_manhattan accumulates chromosome positions", {
  df <- data.frame(
    CHR = c(1, 1, 2, 2, 3, 3),
    BP = c(100, 500, 100, 900, 200, 800),
    P = c(0.5, 0.001, 0.3, 0.01, 0.4, 0.05)
  )
  plot <- plot_manhattan(df, title = "Manhattan")
  expect_s3_class(plot, "ggplot")
  expect_equal(ggplot2::get_labs(plot)$title, "Manhattan")
  tmp <- tempfile(fileext = ".pdf")
  plot_manhattan(df, output_file = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_qq computes lambda and draws the identity line", {
  df <- data.frame(P = runif(200, 0.001, 1))
  plot <- plot_qq(df, p_column = "P")
  expect_s3_class(plot, "ggplot")
  expect_match(ggplot2::get_labs(plot)$title, "lambda")
  expect_true(any(vapply(plot$layers,
                         function(l) inherits(l$geom, "GeomAbline"), logical(1))))
})

test_that("plot_manhattan and plot_qq validate inputs", {
  df <- data.frame(CHR = 1, BP = 100, P = 0.5)
  expect_error(plot_manhattan(df, chr_column = "nope"), "not a column")
  expect_error(plot_manhattan(data.frame(CHR = "x", BP = "y", P = "z")),
               "must be numeric")
  expect_error(plot_qq(df, p_column = "nope"), "not a column")
  expect_error(plot_qq(data.frame(P = c(0, 0))), "positive")
})
test_that("plot_sankey draws an alluvial/sankey diagram", {
  df <- data.frame(
    source = c("a", "a", "b", "b"),
    target = c("x", "y", "x", "y"),
    value = c(10, 5, 8, 12)
  )
  plot <- plot_sankey(df, title = "Sankey")
  expect_s3_class(plot, "ggplot")
  expect_equal(ggplot2::get_labs(plot)$title, "Sankey")
  tmp <- tempfile(fileext = ".pdf")
  plot_sankey(df, output_file = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_sankey supports curve types and custom fill", {
  df <- data.frame(source = c("a", "b"), target = c("x", "y"),
                   value = c(5, 3), group = c("g1", "g2"))
  plot <- plot_sankey(df, curve_type = "sigmoid", fill_column = "group")
  expect_s3_class(plot, "ggplot")
  expect_error(plot_sankey(df, curve_type = "bogus"), "should be one of")
})

test_that("plot_sankey validates inputs", {
  df <- data.frame(source = "a", target = "x", value = 1)
  expect_error(plot_sankey(df, source_column = "nope"), "not a column")
  expect_error(plot_sankey(data.frame(source = "a", target = "x", value = "z")),
               "must be numeric")
  expect_error(plot_sankey(df, flow_alpha = 2), "between 0 and 1")
})


test_that("plot_circular_heatmap draws a gene-by-sample ring heatmap", {
  matrix <- matrix(
    rnorm(40 * 8, mean = 8, sd = 1.5),
    nrow = 40,
    dimnames = list(paste0("Gene", 1:40), paste0("Sample", 1:8))
  )
  scaled <- plot_circular_heatmap(
    matrix, scale = "row", cluster = TRUE, dend_side = "inside",
    rownames_side = "outside", show_colnames = TRUE,
    draw_plot = FALSE
  )
  expect_equal(dim(scaled), c(40, 8))
  expect_true(all(is.finite(scaled)))
  # row Z-scores: near-zero mean and unit sd per row
  expect_equal(round(unname(rowMeans(scaled)), 6), rep(0, 40))
  expect_equal(round(unname(apply(scaled, 1, stats::sd)), 6), rep(1, 40))
})

test_that("plot_circular_heatmap saves PDF output", {
  matrix <- matrix(
    rnorm(30 * 6),
    nrow = 30,
    dimnames = list(paste0("Gene", 1:30), paste0("Sample", 1:6))
  )
  output <- tempfile(fileext = ".pdf")
  plot_circular_heatmap(
    matrix, scale = "row", title = "Circular test",
    output_file = output, figure_width = 8, figure_height = 8
  )
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 1000)
})

test_that("plot_circular_heatmap validates inputs", {
  matrix <- matrix(rnorm(12), nrow = 4, dimnames = list(paste0("G", 1:4), paste0("S", 1:3)))
  expect_error(plot_circular_heatmap(matrix, dend_side = "sideways"), "should be one of")
  expect_error(plot_circular_heatmap(matrix, rownames_side = "sideways"), "should be one of")
  expect_error(
    plot_circular_heatmap(matrix, dend_side = "inside", rownames_side = "inside"),
    "must differ"
  )
  expect_error(
    plot_circular_heatmap(matrix, heatmap_colors = c("red", "blue")),
    "exactly three"
  )
  expect_error(plot_circular_heatmap(matrix(rep(NA_real_, 4), nrow = 2)), "no finite")
})

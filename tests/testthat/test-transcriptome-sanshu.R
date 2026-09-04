test_that("transcriptome summary plots return ggplot objects", {
  style <- choose_plot_style()
  counts <- data.frame(
    Comparison = rep(c("A_vs_B", "C_vs_B"), 2),
    Direction = rep(c("Up", "Down"), each = 2),
    Count = c(12, 9, 7, 5)
  )
  expect_s3_class(deg_count_bar(counts, style = style), "ggplot")
  expect_s3_class(
    deg_upset(list(A = letters[1:8], B = letters[5:12], C = letters[c(2, 6, 10)]),
               style = style),
    "ggplot"
  )
})

test_that("matrix names are checked instead of silently reordered", {
  correlation <- diag(3)
  rownames(correlation) <- c("S1", "S2", "S3")
  colnames(correlation) <- c("S2", "S1", "S3")
  expect_error(sample_correlation_heatmap(correlation), "same order")
})

test_that("GO DAG rejects cycles", {
  nodes <- data.frame(
    ID = c("GO:1", "GO:2"), Term = c("one", "two"),
    Pvalue = c(0.01, 0.02), GeneRatio = c("2/10", "3/10")
  )
  edges <- data.frame(Parent = c("GO:1", "GO:2"), Child = c("GO:2", "GO:1"))
  expect_error(go_dag(nodes, edges), "no root|cycle")
})

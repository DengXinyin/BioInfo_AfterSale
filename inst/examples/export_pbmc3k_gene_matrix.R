# Export a gene-by-cell normalized expression matrix for Python gene UMAP.
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required.")

suppressPackageStartupMessages(library(Seurat))
set.seed(20260825)

input_rds <- "/Users/hsinyinteng/Multi_Omics/scRNA_seq/Seurat/pbmc3k_final.rds"
matrix_file <- "/private/tmp/pbmc3k_gene_by_cell.mtx"
gene_file <- "/private/tmp/pbmc3k_gene_names.tsv"
cell_file <- "/private/tmp/pbmc3k_cell_groups.tsv"

object <- readRDS(input_rds)
cells_1 <- sample(WhichCells(object, idents = "CD14+ Mono"), 400L)
cells_2 <- sample(WhichCells(object, idents = "Naive CD4 T"), 400L)
selected_cells <- c(cells_1, cells_2)
subset_object <- subset(object, cells = selected_cells)

expression <- GetAssayData(subset_object, assay = "RNA", layer = "data")
expression <- expression[, selected_cells, drop = FALSE]
Matrix::writeMM(expression, matrix_file)
utils::write.table(
  data.frame(gene = rownames(expression)), gene_file,
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  data.frame(
    cell = colnames(expression),
    celltype = as.character(Idents(subset_object)[colnames(expression)])
  ),
  cell_file, sep = "\t", quote = FALSE, row.names = FALSE
)

cat("matrix_rows_genes=", nrow(expression), "\n", sep = "")
cat("matrix_columns_cells=", ncol(expression), "\n", sep = "")
cat("nonzero_values=", length(expression@x), "\n", sep = "")
cat("matrix=", matrix_file, "\n", sep = "")
cat("genes=", gene_file, "\n", sep = "")
cat("cells=", cell_file, "\n", sep = "")

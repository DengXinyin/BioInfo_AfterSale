## Generate the four reference figures used by the plotting tutorials.
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "R", "utils.R"))
source(file.path(root, "R", "choose_plot_style.R"))
source(file.path(root, "R", "PCA.R"))
source(file.path(root, "R", "Volcano.R"))
source(file.path(root, "R", "Venn.R"))
source(file.path(root, "R", "GSEA.R"))
source(file.path(root, "R", "GO_KEGG_plot.R"))
source(file.path(root, "R", "Heatmap.R"))

dir.create(file.path(root, "docs", "images"), recursive = TRUE, showWarnings = FALSE)
style <- choose_plot_style(
  # All figure text prefers Times New Roman.
  font_family = "Times New Roman", theme = "bw", dpi = 400,
  figure_width = 9, figure_height = 7,
  title = list(size = 24, bold = TRUE),
  axis_title = list(size = 18, bold = TRUE),
  axis_text = list(size = 16, bold = TRUE),
  legend_title = list(size = 16, bold = TRUE),
  legend_text = list(size = 15, bold = TRUE),
  data_label = list(size = 14, bold = TRUE),
  legend = list(position = "right", frame = FALSE),
  panel = list(border = TRUE, major_grid = TRUE, minor_grid = FALSE),
  group_palette = c("#335372", "#E25659", "#A3A4CA")
)

set.seed(2026)
pca_data <- data.frame(
  Sample = paste0(rep(c("Control", "TA0", "TA4"), each = 5), "-", 1:5),
  PC1 = c(rnorm(5, 0, 12), rnorm(5, 38, 13), rnorm(5, -34, 15)),
  PC2 = c(rnorm(5, -5, 8), rnorm(5, 2, 9), rnorm(5, 15, 10)),
  Group = rep(c("Control", "TA0%", "TA4%"), each = 5)
)
plot_pca(pca_data, pc1 = "PC1", pc2 = "PC2", group = "Group", sample = "Sample",
         variance = c(77.87, 7.60), ellipse = TRUE, show_labels = TRUE,
         legend_inside = c(0.98, 0.98),
         title = "Principal Component Analysis", style = style,
         output_file = file.path(root, "docs", "images", "PCA.png"))

n_gene <- 5000
deg_data <- data.frame(
  GeneID = paste0("Gene", seq_len(n_gene)),
  log2FC = rnorm(n_gene, 0, 1.8),
  pvalue = pmin(runif(n_gene), 10^(-rexp(n_gene, rate = 0.8)))
)
plot_volcano(deg_data, log2fc = "log2FC", pvalue = "pvalue", fc_cutoff = 1,
             p_cutoff = 0.05, status_colors = c("#E25659", "#D4D4D4", "#335372"),
             title = "TA4% vs Control", legend_title = NULL, style = style,
             output_file = file.path(root, "docs", "images", "Volcano.png"))

gene_universe <- paste0("Gene", 1:5000)
venn_sets <- list(
  `TA4% vs TA0%` = sample(gene_universe, 1800),
  `TA4% vs Control` = sample(gene_universe, 2200),
  `TA0% vs Control` = sample(gene_universe, 1500)
)
plot_venn(venn_sets, colors = c("#A3A4CA", "#E25659", "#335372"),
          title = "Differentially Expressed Gene Overlap", style = style,
          output_file = file.path(root, "docs", "images", "Venn.png"),
          width = 8, height = 7)

n_rank <- 12000
ranked_metric <- sort(rnorm(n_rank), decreasing = TRUE)
hit_positions <- sort(sample(seq_len(n_rank), 260))
is_hit <- seq_len(n_rank) %in% hit_positions
hit_weight <- abs(ranked_metric) * is_hit
running_score <- cumsum(hit_weight / sum(hit_weight) - (!is_hit) / sum(!is_hit))
plot_gsea(running_score = running_score, hits = hit_positions, ranked_metric = ranked_metric,
          title = "Peptide cross-linking",
          statistics = c(NES = 2.22, `P value` = 6.1e-09, `Adjusted P` = 2.8e-06),
          colors = c("#E25659", "#335372"), style = style,
          output_file = file.path(root, "docs", "images", "GSEA.png"),
          width = 9.2, height = 6.8)

enrichment_table <- data.frame(
  Description = c("DNA repair", "Cell cycle checkpoint", "Response to oxidative stress",
                  "Mitochondrial organization", "Protein folding", "RNA processing",
                  "Apoptotic signaling pathway", "Chromatin organization", "Lipid metabolic process",
                  "Vesicle-mediated transport"),
  pvalue = c(.00008, .0002, .0007, .0015, .003, .006, .009, .014, .021, .032),
  p.adjust = c(.001, .0018, .0045, .008, .012, .019, .026, .034, .041, .048),
  RichFactor = c(.42, .38, .51, .33, .47, .29, .44, .36, .31, .27),
  Count = c(18, 16, 14, 13, 12, 11, 10, 9, 9, 8), stringsAsFactors = FALSE
)
plot_args <- list(result = enrichment_table, plot_type = "dotplot", filter_by = "p.adjust",
                  cutoff = .05, show_category = 10, label = "Description", size = "Count",
                  style = style, figure_width = 9, figure_height = 6.5, dpi = 300)
do.call(GO_KEGG_plot, c(plot_args, list(x = "pvalue", x_transform = "neg_log10",
  color = "RichFactor", order_by = "pvalue", decreasing = TRUE,
  x_label = expression(-log[10](Pvalue)), color_label = "Rich factor",
  size_label = "Gene count", output_file = file.path(root, "docs", "images", "GO_KEGG_plot_pvalue.png"))))
do.call(GO_KEGG_plot, c(plot_args, list(x = "p.adjust", x_transform = "neg_log10",
  color = "RichFactor", order_by = "p.adjust", decreasing = TRUE,
  x_label = expression(-log[10](Padj)), color_label = "Rich factor",
  size_label = "Gene count", output_file = file.path(root, "docs", "images", "GO_KEGG_plot_padj.png"))))
do.call(GO_KEGG_plot, c(plot_args, list(x = "RichFactor", color = "p.adjust",
  color_transform = "neg_log10", order_by = "RichFactor", decreasing = TRUE,
  x_label = "Rich factor", color_label = expression(-log[10](Padj)), size_label = "Gene count",
  output_file = file.path(root, "docs", "images", "GO_KEGG_plot_richfactor.png"))))

set.seed(20260816)
sample_ids <- c(paste0("Control_", 1:4), paste0("Treatment_", 1:4))
gene_ids <- c(paste0("Up_", sprintf("%02d", 1:8)), paste0("Down_", sprintf("%02d", 1:8)),
              paste0("Stable_", sprintf("%02d", 1:4)))
expression_matrix <- matrix(rnorm(length(gene_ids) * length(sample_ids), 8, .45),
  nrow = length(gene_ids), dimnames = list(gene_ids, sample_ids))
expression_matrix[1:8, 5:8] <- expression_matrix[1:8, 5:8] + 2.2
expression_matrix[9:16, 5:8] <- expression_matrix[9:16, 5:8] - 2.2
sample_group <- setNames(rep(c("Control", "Treatment"), each = 4), sample_ids)
Heatmap_plot(expression_matrix, group = sample_group,
  group_colors = c(Control = "#4E79A7", Treatment = "#E15759"), scale = "row",
  cluster_rows = TRUE, cluster_columns = FALSE, column_names_rot = 45,
  row_names_font_family = "Times New Roman",
  column_names_font_family = "Times New Roman", title = "Synthetic treatment-response genes",
  title_font_family = "Times New Roman", output_file = file.path(root, "docs", "images", "Heatmap_plot.png"),
  figure_width = 8, figure_height = 7, dpi = 300)

message("Reference figures written to ", file.path(root, "docs", "images"))

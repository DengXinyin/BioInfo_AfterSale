# BioInfoAfterSale 0.3.0

- Added `plot_pca()`: principal-component analysis score plot with per-group
  confidence ellipses, sample labels, configurable PC pair, and variance-explain
  axis labels.
- Added `plot_venn()`: two/three-set Venn diagram from a 0/1 indicator data
  frame or a named list.
- Added `plot_violin_box()`: violin/box/jitter distribution plot with optional
  faceting.
- Added `plot_scatter()`: scatter plot with trend line, correlation annotation,
  grouping, faceting, and optional marginal density/box layers.
- Added `plot_distribution()`: histogram, density, or overlaid distribution plot
  with optional grouping.
- Added `plot_survival()`: Kaplan-Meier curve with confidence intervals,
  log-rank P value, risk table, and optional covariate faceting.
- Added `plot_forest()`: Cox regression forest plot from a `coxph` model or a
  precomputed hazard-ratio table.
- Added `plot_manhattan()`: GWAS Manhattan plot with accumulated chromosome
  positions and Bonferroni/suggestive threshold annotations.
- Added `plot_qq()`: Q-Q plot of observed versus expected -log10 P values with
  genomic inflation factor (lambda).
- Added `plot_sankey()`: alluvial/Sankey diagram from a long-form edge table.
- All `plot_*` functions accept a `style` object from `choose_plot_style()` and
  optional PDF/PNG `output_file` output, following the package conventions.

# BioInfoAfterSale 0.2.0

- Added configurable X/Y axis labels, X-axis limits, breaks, and expansion.
- Added point-size range, point opacity, and explicit color/size guide order.
- Added reusable panel border, major/minor grid, and legend-box style controls.
- Changed `size_breaks = NULL` to show at most three representative observed
  Count values in dot-plot legends; explicit `size_breaks` remain supported.
- Added `highlight_terms`, `highlight_color`, and `highlight_bold` to emphasize
  selected pathway labels in both data-frame and native `enrichResult` plots.
- Fixed the default color mapping for unfiltered native enrichment results so
  it uses `p.adjust` instead of treating `GeneRatio` as a color variable.
- Retained selectable output width and height in inches.
- Added `Heatmap_plot()` based on ComplexHeatmap, with configurable labels,
  clustering, Z-score scaling, sample groups, ordered legends, title styling,
  italic gene labels, column-label angles, and PDF/PNG output.
# BioInfoAfterSale 0.2.1

* Added independent LEfSe analysis, LDA score and cladogram functions in
  `R/LEFSE.R`, extracted from the `micro:v2.45` microbiomeMarker workflow.
* Added a modular WGCNA workflow derived from
  `192.168.30.202:23099/wgcna/wgcna:v1.6.7`, including expression QC,
  soft-power selection, module construction, no-trait module-sample summaries,
  optional module-trait correlations, kME, TOM network export, result writing,
  and diagnostic plotting. All public entry points use the `WGCNA_` prefix.
* Documented per-function R package requirements and added an executable
  WGCNA tutorial using the validated 935-gene by 15-sample dataset.

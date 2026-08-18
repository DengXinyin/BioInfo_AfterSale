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

* Added a modular WGCNA workflow derived from
  `192.168.30.202:23099/wgcna/wgcna:v1.6.7`, including expression QC,
  soft-power selection, module construction, no-trait module-sample summaries,
  optional module-trait correlations, kME, TOM network export, result writing,
  and diagnostic plotting. All public entry points use the `WGCNA_` prefix.
* Documented per-function R package requirements and added an executable
  WGCNA tutorial using the validated 935-gene by 15-sample dataset.

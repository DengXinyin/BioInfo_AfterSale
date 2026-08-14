# BioInfoAfterSale 0.2.0

- Added configurable X/Y axis labels, X-axis limits, breaks, and expansion.
- Added point-size range, point opacity, and explicit color/size guide order.
- Added reusable panel border, major/minor grid, and legend-box style controls.
- Fixed `size_breaks = NULL` so Count legends use automatic ggplot2 breaks.
- Retained selectable output width and height in inches.
- Added `Heatmap_plot()` based on ComplexHeatmap, with configurable labels,
  clustering, Z-score scaling, sample groups, ordered legends, title styling,
  italic gene labels, column-label angles, and PDF/PNG output.

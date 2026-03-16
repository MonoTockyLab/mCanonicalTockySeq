# Plot Tocky-Time Heatmap

Generates a binned and smoothed heatmap ordered by Tocky Time.

## Usage

``` r
mPlotTockyHeatmap(object, seurat_obj, genes, ...)

# S4 method for class 'mCanonicalTockyObj'
mPlotTockyHeatmap(
  object,
  seurat_obj,
  genes,
  n_bins = 100,
  ordering_method = c("peak", "cluster"),
  span = 0.5,
  scale_rows = TRUE
)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- seurat_obj:

  Original Seurat object.

- genes:

  Character vector of genes to plot.

- ...:

  Additional arguments passed to methods.

- n_bins:

  Integer. Number of bins (default 100).

- ordering_method:

  Character. "peak" (sort by max time) or "cluster".

- span:

  Numeric. The smoothing parameter (alpha) for LOESS smoothing across
  Tocky Time bins. Default is 0.5.

- scale_rows:

  Logical. Whether to perform row-wise Z-score scaling on the smoothed
  expression matrix. Default is TRUE.

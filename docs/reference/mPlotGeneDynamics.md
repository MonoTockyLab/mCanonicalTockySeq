# Plot Gene Expression Dynamics by Group

Visualizes single-cell gene expression trends along the Tocky Time axis.
Accepts multiple genes and automatically generates a grid layout.
Automatically detects Lineage_Bias from the mCanonicalTockyObj.

## Usage

``` r
mPlotGeneDynamics(object, seurat_obj, features, ...)

# S4 method for class 'mCanonicalTockyObj'
mPlotGeneDynamics(
  object,
  seurat_obj,
  features,
  group_by = "Lineage_Bias",
  span = 0.8,
  jitter_amount = 0.2,
  pt_alpha = 0.2,
  m = 1.15,
  ncol = 2
)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- seurat_obj:

  Original Seurat object.

- features:

  Character vector of genes to plot.

- ...:

  Additional arguments passed to methods.

- group_by:

  Optional metadata column in Seurat or TockyData to color by. Defaults
  to "Lineage_Bias".

- span:

  Numeric. Smoothing span for LOESS (default 0.8).

- jitter_amount:

  Numeric. Amount of vertical jitter applied to expression values to
  reduce overplotting (default 0.2).

- pt_alpha:

  Numeric. Transparency level for the plotted cells, ranging from 0
  (transparent) to 1 (opaque) (default 0.2).

- m:

  Numeric. Multiplier for the maximum y-axis limit to ensure curves and
  data points fit comfortably within the plot (default 1.15).

- ncol:

  Integer. Number of columns for the plot grid (default 2).

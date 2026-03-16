# Launch Interactive 3D Fate Explorer App

Launches a local Shiny web application to explore the 3D Tocky
trajectory. Utilizes dynamically anchored overlapping multi-identity
tubes from mExtractTrajectoryTubes().

## Usage

``` r
LaunchTocky3DApp(
  object,
  seurat_obj,
  window_size = 10,
  step_size = 1,
  span = 0.4
)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- seurat_obj:

  Original Seurat object to fetch gene expression from.

- window_size:

  Numeric. Width of the sliding angle slice (default 10).

- step_size:

  Numeric. Step size for the sliding window (default 1).

- span:

  Numeric. LOESS smoothing span (default 0.4).

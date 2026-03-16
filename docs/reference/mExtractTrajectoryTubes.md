# Extract Independent Biological Trajectory Tubes

Models independent developmental trajectories for an arbitrary number of
genes. For each gene, it isolates the "Core" expressing cells to
calculate an expression-weighted barycenter and covariance matrix,
bypassing the gravitational pull of ambient RNA noise. It then draws a
perfect geometric ellipse using the Chi-Square distribution. Cells
falling within the confidence interval of a gene's trajectory are marked
as TRUE for that specific tube, allowing for overlapping dual
identities.

## Usage

``` r
mExtractTrajectoryTubes(object, ...)

# S4 method for class 'mCanonicalTockyObj'
mExtractTrajectoryTubes(
  object,
  seurat_obj,
  features,
  ci = 0.9,
  window_size = 10,
  span = 0.4,
  reset = TRUE
)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- ...:

  Additional arguments passed to methods.

- seurat_obj:

  Original Seurat object.

- features:

  Character vector of genes to model trajectories for (e.g., c("Cd4",
  "Cd8a")).

- ci:

  Numeric. Confidence interval for the ellipses (default 0.90).

- window_size:

  Numeric. Width of the sliding angle slice (default 10).

- span:

  Numeric. LOESS smoothing span for the barycenters (default 0.4).

- reset:

  Logical. If TRUE (default), deletes previously created lineage tubes.
  If FALSE, appends new tubes and updates the composite
  Lineage_Identity.

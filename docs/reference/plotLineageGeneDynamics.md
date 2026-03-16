# Plot Gene Dynamics Trends Along Trajectory Tubes

Generates a high-clarity plot showing smoothed maturation trends for
specific genes. Automatically maps cells to their respective overlapping
biological tubes (e.g., Tube_Cd4, Tube_Cd8a) and plots the expression
dynamics across Tocky Time. Cells with dual identities properly
contribute to the shared origins of multiple trajectories.

## Usage

``` r
plotLineageGeneDynamics(object, seurat_obj, features, ...)

# S4 method for class 'mCanonicalTockyObj'
plotLineageGeneDynamics(object, seurat_obj, features, ncol = 2)
```

## Arguments

- object:

  An mCanonicalTockyObj containing tube assignments.

- seurat_obj:

  The original Seurat object to fetch gene expression from.

- features:

  Character vector of genes to plot.

- ...:

  Additional arguments passed to methods.

- ncol:

  Integer. Number of columns for the plot grid (default 2).

# Plot Canonical Tocky Locus Features

Generates a matrix plot showing Tocky progression and gene expression
across binned maturation stages in fate space.

## Usage

``` r
plotTockyLocusFeatures(object, seurat_obj, features, ...)

# S4 method for class 'mCanonicalTockyObj'
plotTockyLocusFeatures(object, seurat_obj, features, bins = 5)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- seurat_obj:

  The original Seurat object to fetch expression from.

- features:

  Character vector of genes to plot.

- ...:

  Additional arguments passed to methods.

- bins:

  Number of angle bins (default 5).

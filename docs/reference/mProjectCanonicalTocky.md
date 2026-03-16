# Project New ScRNA-seq Data into an Established Canonical Tocky Space

Performs Redundancy Analysis (RDA) on a new expression matrix (e.g.,
human orthologs) using the exact explanatory variable matrix (Z) and
internal mappings from a reference Tocky model (e.g., mouse).

## Usage

``` r
mProjectCanonicalTocky(X_mat, ref_obj, seurat_obj, ...)

# S4 method for class 'Matrix'
mProjectCanonicalTocky(X_mat, ref_obj, seurat_obj, scale_Z = TRUE)
```

## Arguments

- X_mat:

  A Matrix (dense or sparse) of expression data for the new dataset.

- ref_obj:

  The reference mCanonicalTockyObj.

- seurat_obj:

  The Seurat object corresponding to the cells in X_mat.

- ...:

  Additional arguments passed to methods.

- scale_Z:

  Logical. Whether to scale the Z matrix (default TRUE).

## Value

An mCanonicalTockyObj containing the projected cells.

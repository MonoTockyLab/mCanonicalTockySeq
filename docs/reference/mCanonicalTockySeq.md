# Canonical Redundancy Analysis for Tocky Differentiation

Canonical Redundancy Analysis for Tocky Differentiation

## Usage

``` r
mCanonicalTockySeq(
  object,
  temporal_col,
  b_ident = "B",
  br_ident = "BR",
  r_ident = "R",
  lineage_col,
  lineage_idents,
  lineage_names = NULL,
  top_n = 100,
  custom_genes = NULL,
  layer = "data",
  scale_Z = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- temporal_col:

  Character. Metadata column for temporal anchors.

- b_ident:

  Character. Identity for the Blue (New) landmark (default "B").

- br_ident:

  Character. Identity for the Blue-Red (Persistent) landmark (default
  "BR").

- r_ident:

  Character. Identity for the Red (Arrested) landmark (default "R").

- lineage_col:

  Character. Metadata column for lineage endpoints.

- lineage_idents:

  Character vector of identities in Seurat object for lineage endpoints.

- lineage_names:

  Optional character vector. Biological names for the lineages (e.g.,
  c("CD4", "CD8")).

- top_n:

  Integer. Number of marker genes to extract per group (default 100).

- custom_genes:

  Optional character vector to skip marker calculation.

- layer:

  Character. The Seurat layer to use (default "data").

- scale_Z:

  Logical. Whether to scale the Z matrix (default TRUE).

## Value

An object of class `mCanonicalTockyObj`.

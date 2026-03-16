# S4 Class to hold mCanonicalTockySeq Results

S4 Class to hold mCanonicalTockySeq Results

## Slots

- `X`:

  Matrix. The single-cell expression matrix (Genes x Cells), accepts
  either dense or sparse matrices.

- `Z`:

  matrix. The reference landmark matrix (Genes x Vectors).

- `metadata`:

  data.frame. Cell-level metadata.

- `genes`:

  character. The features used for the RDA.

- `expression_scores`:

  data.frame. Gene loadings.

- `fitted_cell_scores`:

  data.frame. Fitted cell scores.

- `cell_scores`:

  data.frame. Projected cell scores.

- `biplot`:

  data.frame. Constraint vectors.

- `marker_info`:

  list. Provenance of landmarks.

- `TockyData`:

  data.frame. Gradient mapping and lineage trajectory results.

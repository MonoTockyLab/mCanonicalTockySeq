# Run TradeSeq Differential Dynamics Analysis (Optimized for 2 Lineages)

Evaluates differential gene expression across exactly two developmental
tubes. Employs aggressive proportional thresholding to drop
uninformative genes, bypassing severe RAM and processing bottlenecks in
tradeSeq.

## Usage

``` r
mRunTradeSeqTocky(object, seurat_obj, ...)

# S4 method for class 'mCanonicalTockyObj'
mRunTradeSeqTocky(
  object,
  seurat_obj,
  tubes = NULL,
  min_pct = 0.05,
  n_knots = 5,
  n_cores = 1,
  p_adjust_method = "fdr"
)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- seurat_obj:

  The Seurat object to fetch counts from.

- ...:

  Additional arguments passed to methods.

- tubes:

  Character vector of exactly two genes defining the tubes to compare
  (e.g., c("Cd4", "Cd8a")). If NULL, auto-detects if exactly two exist.

- min_pct:

  Numeric. Minimum fraction of cells a gene must be expressed in
  (default 0.05, i.e., 5%).

- n_knots:

  Integer. Number of knots for the GAM (default 5).

- n_cores:

  Integer. Number of cores for parallel processing (default 1).

- p_adjust_method:

  Character. Correction method: "fdr" (default), "bonferroni", "holm",
  etc.

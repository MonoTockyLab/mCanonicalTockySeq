# Plot Canonical Tocky Ordination with Continuous Gradient Legend

Visually evaluates the constrained ordination space. Automatically
detects the number of available RDA axes and generates a multi-panel
layout showing all sequential axis pairs (1v2, 2v3, 3v4, etc.) followed
by a continuous Tocky Time legend.

## Usage

``` r
mPlotCanonicalTocky(object, ...)

# S4 method for class 'mCanonicalTockyObj'
mPlotCanonicalTocky(object, alpha_level = 0.2)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- ...:

  Additional arguments passed to methods.

- alpha_level:

  Numeric transparency level for cells (0-1). Default is 0.2.

# Gradient Tocky Sequence (Slerp Model)

Projects cells onto a piecewise Slerp manifold defined by the landmarks
B -\> BR -\> R. Automatically extracts the landmark vectors from the
provided `mCanonicalTockyObj`. Cells with negative correlations to all
three landmarks are classified as "Timer Negative" (NA), as they lie in
the ordination space opposite to the Tocky trajectory.

## Usage

``` r
mGradientTockySeq(object, ...)

# S4 method for class 'mCanonicalTockyObj'
mGradientTockySeq(object, filter_negative = TRUE, intensity_q = 0.75)
```

## Arguments

- object:

  An object of class `mCanonicalTockyObj`.

- ...:

  Additional arguments.

- filter_negative:

  Logical. If TRUE, assigns NA to cells with negative projection to all
  three landmarks.

- intensity_q:

  Numeric. Quantile threshold (0 to 1) for normalized intensity. Cells
  falling below this threshold are treated as non-participating and
  assigned NA for angle. Default is 0.75.

## Value

An updated `mCanonicalTockyObj` with the Slerp mapping results stored in
the `@TockyData` slot.

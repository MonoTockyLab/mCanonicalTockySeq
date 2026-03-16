# Calculate Fate Commitment Scores

Calculates the dot product projection of cells towards the lineage
endpoints. Automatically detects the lineage destinations established
during mCanonicalTockySeq().

## Usage

``` r
mGetFateScores(object, ...)

# S4 method for class 'mCanonicalTockyObj'
mGetFateScores(object)
```

## Arguments

- object:

  An mCanonicalTockyObj.

- ...:

  Additional arguments passed to methods.

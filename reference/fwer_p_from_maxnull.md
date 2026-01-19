# Compute adjusted p-values from max null distribution

Helper function to compute FWER-corrected p-values given a max null
distribution.

## Usage

``` r
fwer_p_from_maxnull(observed, max_null)
```

## Arguments

- observed:

  Numeric vector of observed values

- max_null:

  Numeric vector of max statistics under null

## Value

Numeric vector of adjusted p-values

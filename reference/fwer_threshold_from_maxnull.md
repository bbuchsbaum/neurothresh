# Compute discrete permutation cutoff from max null distribution

Returns a cutoff `u` such that the standard plus-one corrected
permutation rule rejects for `observed > u` at level `alpha`.

## Usage

``` r
fwer_threshold_from_maxnull(max_null, alpha)
```

## Arguments

- max_null:

  Numeric vector of max statistics under null

- alpha:

  Significance level in (0, 1)

## Value

Numeric scalar cutoff

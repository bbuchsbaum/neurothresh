# Empirical FDR for LISA filtered maps

Estimates local FDR/q-values in the filtered-statistic space using the
empirical filtered null distribution and the empirical observed-map
mixture distribution.

## Usage

``` r
lisa_fdr_map(
  filtered_obs,
  filtered_nulls,
  mask = NULL,
  tail = c("pos", "neg", "two"),
  q = 0.05,
  p0 = 1,
  monotone = TRUE
)
```

## Arguments

- filtered_obs:

  Numeric observed filtered statistic vector or volume.

- filtered_nulls:

  Null filtered statistics. Accepts a matrix with one row per
  permutation and one column per voxel, or a list of numeric vectors.

- mask:

  Optional logical mask when `filtered_obs` is a volume.

- tail:

  One of `"pos"`, `"neg"`, or `"two"`.

- q:

  Target FDR level used to compute `threshold` and `reject`.

- p0:

  Conservative null proportion. Defaults to 1.

- monotone:

  Whether to enforce monotone q-values over evidence strength.

## Value

List with vector `fdr`, scalar `threshold`, logical `reject`, and an
empirical `curve` data frame.

# Canonicalize a statistic map to a Z-equivalent volume

Convenience wrapper around
[`canonicalize_stat()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize_stat.md)
that returns only the Z-equivalent statistic map. This is the canonical
first step before calling inference functions in neurothresh.

## Usage

``` r
canonicalize(
  stat_vol,
  stat_type = c("Z", "t", "neglog10p"),
  df = NULL,
  tail = c("pos", "neg", "two"),
  p_side = c("one", "two"),
  sign_vol = NULL
)
```

## Arguments

- stat_vol:

  A NeuroVol object containing the statistical map

- stat_type:

  Character string specifying the input type:

  "Z"

  :   Z-scores (identity transform)

  "t"

  :   t-statistics (requires `df`)

  "neglog10p"

  :   Negative log10 p-values

- df:

  Degrees of freedom, required when `stat_type = "t"`

- tail:

  Character string specifying the tail:

  "pos"

  :   Positive tail (activations)

  "neg"

  :   Negative tail (deactivations)

  "two"

  :   Two-sided test

- p_side:

  Character string for -log10(p) inputs specifying whether the p-values
  are one-sided or two-sided: "one" or "two"

- sign_vol:

  Optional NeuroVol providing sign information for -log10(p) inputs to
  recover directionality

## Value

The Z-equivalent statistic volume (same type as `stat_vol`). For richer
output (one-sided/two-sided p-values, mask), call
[`canonicalize_stat()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize_stat.md)
directly.

## Examples

``` r
if (FALSE) { # \dontrun{
zvol <- canonicalize(t_map, stat_type = "t", df = 24)
res <- hier_scan(zvol, parcels = atlas, alpha = 0.05)
} # }
```

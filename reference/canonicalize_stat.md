# Canonicalize statistical maps to Z-scores

Converts various statistical map types (Z-scores, t-statistics, or
-log10(p) values) to a common Z-score representation for unified
thresholding procedures.

## Usage

``` r
canonicalize_stat(
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

A list with components:

- Zeq:

  NeuroVol of equivalent Z-scores

- p_one:

  NeuroVol of one-sided p-values

- p_two:

  NeuroVol of two-sided p-values

- mask:

  Logical mask of valid voxels

## Details

The conversions are:

- Z-scores: Used directly

- t-statistics: Converted via `qnorm(pt(t, df))`

- -log10(p): Converted via `qnorm(1 - 10^(-x))` with optional sign
  recovery from `sign_vol`

## Examples

``` r
if (FALSE) { # \dontrun{
# From t-statistics with 50 degrees of freedom
result <- canonicalize_stat(t_map, stat_type = "t", df = 50)

# From -log10(p) with sign information
result <- canonicalize_stat(logp_map, stat_type = "neglog10p",
                            sign_vol = t_map)
} # }
```

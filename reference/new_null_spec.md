# Construct a neurothresh null specification

Internal constructor. Factories should call this rather than building
the structure by hand so validation runs uniformly.

## Usage

``` r
new_null_spec(
  stat_vol,
  mask,
  draw,
  n_perm,
  df = NA_real_,
  kind = "unknown",
  meta = list()
)
```

## Arguments

- stat_vol:

  Observed (canonicalized) statistic map.

- mask:

  Logical mask.

- draw:

  `function(b)` returning length-`sum(mask)` numeric.

- n_perm:

  Integer maximum number of draws.

- df:

  Degrees of freedom (numeric scalar, may be `NA`).

- kind:

  Character tag identifying the resampling design.

- meta:

  Optional named list of extras.

## Value

Object of class `"neurothresh_null_spec"`.

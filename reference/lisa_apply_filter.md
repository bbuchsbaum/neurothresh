# Apply a LISA filter

Applies a LISA filter specification to a 3D statistic map. This function
is mostly intended for diagnostics and tests;
[`lisa()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa.md)
applies it internally to observed and null maps.

## Usage

``` r
lisa_apply_filter(
  stat_vol,
  mask = NULL,
  filter = lisa_bilateral(),
  range_scale = NULL
)
```

## Arguments

- stat_vol:

  3D statistic map.

- mask:

  Optional logical 3D mask.

- filter:

  A LISA filter specification from
  [`lisa_bilateral()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_bilateral.md).

- range_scale:

  Optional fixed range scale. If `NULL`, the filter uses the map's
  masked standard deviation.

## Value

Filtered volume-like object aligned to `stat_vol`.

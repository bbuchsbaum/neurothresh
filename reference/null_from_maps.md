# External observed/null statistic maps

Wraps an observed statistic map and externally generated null maps in a
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).
This is the generic bridge for MVPA, searchlight, PLS, ASCA, or any
workflow that already has matched observed and permuted statistic maps.

## Usage

``` r
null_from_maps(
  stat_vol,
  null_maps,
  mask = NULL,
  kind = "external_maps",
  meta = list()
)
```

## Arguments

- stat_vol:

  Observed canonicalized Z-equivalent statistic map.

- null_maps:

  Null maps. Accepts a list of 3D maps, a 4D array whose fourth
  dimension indexes permutations, or a numeric matrix with one row per
  permutation and either `sum(mask)` or `prod(dim(stat_vol))` columns.

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- kind:

  Character kind tag for the resulting
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

- meta:

  Optional named metadata list.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
with deterministic `draw(b)` over supplied maps.

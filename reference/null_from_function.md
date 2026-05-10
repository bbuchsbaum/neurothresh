# Wrap a bare draw function in a null_spec

Convenience for ad-hoc resampling models that already produce a
deterministic-by-index draw function. Most users should reach for the
high-level factories (e.g.
[`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md))
instead.

## Usage

``` r
null_from_function(draw, stat_vol, n_perm, mask = NULL, kind = "custom")
```

## Arguments

- draw:

  `function(b)` returning a length-`sum(mask)` numeric.

- stat_vol:

  Observed statistic map.

- n_perm:

  Number of permutations the draw function supports.

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- kind:

  Optional kind tag (default `"custom"`).

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

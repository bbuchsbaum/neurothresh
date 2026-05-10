# LISA local adaptive inference

Applies post-statistic bilateral filtering to an observed statistic map
and matched null maps, then estimates voxelwise FDR in
filtered-statistic space. This is a local-adaptive companion to
scale-space energy inference.

## Usage

``` r
lisa(
  stat_vol = NULL,
  mask = NULL,
  q = 0.05,
  tail = c("pos", "neg", "two"),
  null_spec = NULL,
  n_perm = NULL,
  seed = NULL,
  filter = lisa_bilateral(),
  p0 = 1,
  monotone = TRUE,
  scale_n_null = 30L,
  backend = c("neuroim2", "reference", "adjoin_sparse"),
  store_filtered_nulls = FALSE
)
```

## Arguments

- stat_vol:

  Optional 3D Z-equivalent statistic map. Defaults to
  `null_spec$stat_vol`. For convenience, `stat_vol` may also be a
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  object, in which case `null_spec` is inferred.

- mask:

  Optional mask override. Defaults to `null_spec$mask`.

- q:

  FDR level.

- tail:

  One of `"pos"`, `"neg"`, or `"two"`.

- null_spec:

  A
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  object. Its `stat_vol`, `mask`, `draw`, and `n_perm` slots define the
  default observed map and null maps.

- n_perm:

  Optional number of null draws to use. Defaults to all draws in
  `null_spec`; may not exceed `null_spec$n_perm`.

- seed:

  Optional random seed set before consuming null draws. Most `null_spec`
  factories are deterministic by draw index, but this keeps LISA aligned
  with other permutation-facing neurothresh methods.

- filter:

  LISA filter specification from
  [`lisa_bilateral()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa_bilateral.md).

- p0:

  Conservative null proportion for the empirical FDR estimator.

- monotone:

  Whether to enforce monotone q-values over evidence strength.

- scale_n_null:

  Number of early null maps used for fixed range-scale estimation when
  `filter$scale = "sd_first_30_nulls"`.

- backend:

  Filtering backend.

- store_filtered_nulls:

  Whether to store the full permutation-by-voxel filtered null matrix in
  `extras`. Defaults to `FALSE`; inference uses streaming tail counts
  either way.

## Value

A `neurothresh_result` with LISA q-values in `pvals` and the filtered
statistic map in `extras$filtered_stat`.

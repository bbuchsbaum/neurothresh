# Scale-Space Localized Energy Inference

Scans a statistic map across a prespecified Gaussian scale-space and
controls FWER over scale, location, and statistic channel using
synchronized null maps.

## Usage

``` r
scale_space_energy(
  stat_vol = NULL,
  mask = NULL,
  scales = scale_space(c(6, 10, 14, 18, 24, 32, 48)),
  statistics = list(local_energy(), signed_coherence(tail = "two")),
  null_fun = NULL,
  null_spec = NULL,
  n_perm = NULL,
  studentize = TRUE,
  alpha = 0.05,
  seed = NULL,
  backend = c("grid_gaussian", "graph_heat")
)
```

## Arguments

- stat_vol:

  3D statistic map, or a
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  when `null_spec` is not supplied separately.

- mask:

  Optional logical analysis mask.

- scales:

  A `scale_space_spec` object.

- statistics:

  List of channel objects such as
  [`local_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/local_energy.md)
  and
  [`signed_coherence()`](https://bbuchsbaum.github.io/neurothresh/reference/signed_coherence.md).

- null_fun:

  Function `function(b)` returning a null statistic draw. It should be
  deterministic by index because calibration reuses draw `b` after
  estimating common observed-plus-null moments.

- null_spec:

  Optional
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  supplying `stat_vol`, `mask`, deterministic `draw(b)`, and `n_perm`.

- n_perm:

  Number of permutations/null draws.

- studentize:

  Whether to studentize transformed fields using null mean and standard
  deviation estimated symmetrically from the observed map plus the null
  draws. Defaults to TRUE.

- alpha:

  FWER level used for summaries.

- seed:

  Optional random seed. Only meaningful when `null_fun` draws on the
  global RNG; the bundled `make_null_fun_*` helpers are deterministic by
  index and ignore this seed.

- backend:

  Smoothing backend. `"grid_gaussian"` uses a separable, mask-normalized
  3D Gaussian convolution. `"graph_heat"` uses a weighted 6-neighbor
  mask graph and applies the continuous-time heat operator `expm(-t L)`
  with `t = sigma_mm^2 / 2`.

## Value

An object of class `scale_space_energy_result`.

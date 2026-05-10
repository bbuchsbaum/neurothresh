# LISA bilateral filter specification

Creates a filter specification for LISA local adaptive inference. The
default backend delegates to
[`neuroim2::bilateral_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.html);
a small reference backend is available for tests and exact mask-aware
semantics on small arrays.

## Usage

``` r
lisa_bilateral(
  radius = 2L,
  sigma_space = 2,
  sigma_range = 2,
  n_iter = 2L,
  scale = c("sd_first_30_nulls", "per_map", "none"),
  boundary = c("neuroim2", "mask_normalized", "lisa2018")
)
```

## Arguments

- radius:

  Integer neighborhood half-width in voxels.

- sigma_space:

  Spatial Gaussian sigma in voxel units.

- sigma_range:

  Range/intensity sigma. With the `reference` backend and a supplied
  range scale, the effective range sigma is `sigma_range * range_scale`.

- n_iter:

  Number of repeated filtering iterations.

- scale:

  Range-scale policy. `"sd_first_30_nulls"` estimates one fixed range
  scale from early null maps. `"per_map"` uses each map's own masked
  standard deviation. `"none"` uses a unit range scale.

- boundary:

  Boundary policy. `"neuroim2"` faithfully delegates to
  [`neuroim2::bilateral_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.html).
  `"mask_normalized"` excludes out-of-mask neighbors in the reference
  backend. `"lisa2018"` is reserved for the paper-style median fallback
  and currently errors.

## Value

Object of class `"neurothresh_lisa_filter"`.

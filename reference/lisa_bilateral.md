# LISA bilateral filter specification

Creates a filter specification for LISA local adaptive inference. LISA
uses the
[`neuroim2::bilateral_filter()`](https://bbuchsbaum.github.io/neuroim2/reference/bilateral_filter.html)
implementation, whose neighborhoods are mask-normalized so only in-mask,
in-bounds neighbors contribute to each filtered statistic.

## Usage

``` r
lisa_bilateral(
  radius = 2L,
  sigma_space = 2,
  sigma_range = 2,
  n_iter = 2L,
  scale = c("sd_first_30_nulls", "per_map", "none")
)
```

## Arguments

- radius:

  Integer neighborhood half-width in voxels.

- sigma_space:

  Spatial Gaussian sigma in voxel units.

- sigma_range:

  Range/intensity sigma. With a fixed range scale, the effective range
  sigma is `sigma_range * range_scale`.

- n_iter:

  Number of repeated filtering iterations.

- scale:

  Range-scale policy. `"sd_first_30_nulls"` estimates one fixed range
  scale from early null maps. `"per_map"` uses each map's own masked
  standard deviation. `"none"` uses a unit range scale.

## Value

Object of class `"neurothresh_lisa_filter"`.

# Monte-Carlo FWHM null (Gaussian-smoothed white noise)

Generates null draws by smoothing white noise with a Gaussian kernel
that yields the requested FWHM. Returns a
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
whose `draw(b)` is length `sum(mask)`.

## Usage

``` r
null_mc_fwhm(template, mask, smoothness, n_perm = 5000, seed = NULL)
```

## Arguments

- template:

  3D template carrying spatial metadata (`NeuroVol`, array, or anything
  with [`dim()`](https://rdrr.io/r/base/dim.html) of length 3).

- mask:

  Logical 3D mask.

- smoothness:

  Either a
  [`smoothness_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_fwhm.md)
  object, or a numeric FWHM in voxel units (length 1 or 3) for backwards
  compatibility.

- n_perm:

  Maximum number of indexed draws supported by the resulting spec.

- seed:

  Optional random seed used to precompute per-index RNG seeds; the
  resulting `draw(b)` is deterministic and restores the caller's RNG
  state.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

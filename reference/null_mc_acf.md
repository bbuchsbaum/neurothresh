# Monte-Carlo ACF null (AFNI-style autocorrelation)

Generates null draws by simulating fields with a mixed Gaussian +
monoexponential autocorrelation model. Returns a
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

## Usage

``` r
null_mc_acf(template, mask, smoothness, n_perm = 5000, seed = NULL)
```

## Arguments

- template:

  3D template carrying spatial metadata (`NeuroVol`, array, or anything
  with [`dim()`](https://rdrr.io/r/base/dim.html) of length 3).

- mask:

  Logical 3D mask.

- smoothness:

  A
  [`smoothness_acf()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_acf.md)
  object, a numeric vector `c(a, b, c)` in millimetres, or a legacy list
  with `a`, `fwhm_vox`, and `lambda_vox`.

- n_perm:

  Maximum number of indexed draws supported by the resulting spec.

- seed:

  Optional random seed used to precompute per-index RNG seeds; the
  resulting `draw(b)` is deterministic and restores the caller's RNG
  state.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

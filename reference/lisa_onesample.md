# LISA one-sample wrapper

LISA one-sample wrapper

## Usage

``` r
lisa_onesample(
  images,
  mask = NULL,
  n_perm = 5000,
  seed = NULL,
  z_equiv = TRUE,
  template = NULL,
  ...
)
```

## Arguments

- images:

  Subject-level images. Accepts any input handled by
  [`as_subject_matrix()`](https://bbuchsbaum.github.io/neurothresh/reference/as_subject_matrix.md):
  list of 3D arrays/`NeuroVol`s, 4D array, `NeuroVec`, `SparseNeuroVec`,
  or a subject-by-voxel matrix (with `template`/`mask`).

- mask:

  Optional logical 3D mask.

- n_perm:

  Number of sign-flip draws.

- seed:

  Optional random seed.

- z_equiv:

  If `TRUE` (default), convert the t statistic to a Z equivalent.

- template:

  Optional template (only used for matrix input).

- ...:

  Additional arguments passed to
  [`lisa()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa.md).

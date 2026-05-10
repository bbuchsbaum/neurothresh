# One-sample subject sign-flip null

Builds an observed one-sample t (or Z-equivalent) map and a
deterministic subject-level sign-flip null indexed by `b`. Each subject
map is multiplied by a single \\\pm 1\\ per draw, preserving
within-subject spatial covariance.

## Usage

``` r
null_one_sample(
  images,
  mask = NULL,
  n_perm = 5000,
  seed = NULL,
  z_equiv = TRUE,
  template = NULL
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

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

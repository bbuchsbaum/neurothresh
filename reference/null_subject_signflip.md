# Subject sign-flip null (alias of one-sample)

The "subject signflip" name describes the mechanism. The factory under
the design name is
[`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md);
this is an alias kept for callers who prefer the mechanism-first name.

## Usage

``` r
null_subject_signflip(
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

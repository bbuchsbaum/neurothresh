# Scale-Space Two-Sample Wrapper

Thin raw-image wrapper around
[`scale_space_energy`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md)
using fixed-size two-group label permutations.

## Usage

``` r
scale_space_twosample(
  subject_images,
  group,
  mask = NULL,
  scales = scale_space(c(6, 10, 14, 18, 24, 32, 48)),
  statistics = list(local_energy(), signed_coherence(tail = "two")),
  n_perm = 1000,
  seed = NULL,
  variance = c("welch", "pooled"),
  z_equiv = TRUE,
  template = NULL,
  alpha = 0.05,
  studentize = TRUE,
  ...
)
```

## Arguments

- subject_images:

  List of subject-level 3D images, or a subject-by-voxel matrix with
  `template`/`mask`.

- group:

  Two-level group label.

- mask:

  Optional logical mask.

- scales:

  Scale-space specification.

- statistics:

  List of scale-space channel objects.

- n_perm:

  Number of label permutations.

- seed:

  Optional random seed.

- variance:

  "welch" (default) or "pooled".

- z_equiv:

  If TRUE, return Z-equivalent maps; otherwise return t maps.

- template:

  Optional template for matrix inputs.

- alpha:

  FWER level.

- studentize:

  Whether to studentize scale-space fields.

- ...:

  Additional arguments passed to
  [`scale_space_energy`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md).

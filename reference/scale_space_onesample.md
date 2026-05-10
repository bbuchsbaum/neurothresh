# Scale-Space One-Sample Wrapper

Thin raw-image wrapper around
[`scale_space_energy`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md)
using subject-level sign-flip resampling.

## Usage

``` r
scale_space_onesample(
  subject_images,
  mask = NULL,
  scales = scale_space(c(6, 10, 14, 18, 24, 32, 48)),
  statistics = list(local_energy(), signed_coherence(tail = "two")),
  n_perm = 1000,
  seed = NULL,
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

- mask:

  Optional logical mask.

- scales:

  Scale-space specification.

- statistics:

  List of scale-space channel objects.

- n_perm:

  Number of sign-flip draws.

- seed:

  Optional random seed.

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

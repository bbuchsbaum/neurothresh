# Subject-Level One-Sample Sign-Flip Null Generator

Builds an observed one-sample statistic map and deterministic
subject-level sign-flip null draws for scale-space and other resampling
methods.

## Usage

``` r
make_null_fun_onesample(
  subject_images,
  mask = NULL,
  n_perm = 1000,
  seed = NULL,
  z_equiv = TRUE,
  template = NULL
)
```

## Arguments

- subject_images:

  List of subject-level 3D images, or a subject-by-voxel matrix with
  `template`/`mask`.

- mask:

  Optional logical mask.

- n_perm:

  Number of sign-flip draws.

- seed:

  Optional random seed.

- z_equiv:

  If TRUE, return Z-equivalent maps; otherwise return t maps.

- template:

  Optional template for matrix inputs.

## Value

List with `stat_vol`, `mask`, deterministic `null_fun`, and metadata.

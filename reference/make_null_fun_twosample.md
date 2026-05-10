# Two-Group Label-Permutation Null Generator

Builds an observed two-group statistic map and deterministic fixed-size
label-permutation null draws.

## Usage

``` r
make_null_fun_twosample(
  subject_images,
  group,
  mask = NULL,
  n_perm = 1000,
  seed = NULL,
  variance = c("welch", "pooled"),
  z_equiv = TRUE,
  template = NULL
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

## Value

List with `stat_vol`, `mask`, deterministic `null_fun`, and metadata.

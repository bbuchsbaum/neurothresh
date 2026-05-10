# Null specification S3 class

A `neurothresh_null_spec` is the uniform handle every permutation-based
inference method consumes. It bundles the observed (canonicalized)
statistic map, the mask used to build the null, a deterministic-by-index
draw function, and metadata describing the resampling design.

## Details

Factories that return `neurothresh_null_spec` objects:
[`null_voxel_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_voxel_signflip.md),
[`null_subject_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_subject_signflip.md),
[`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md),
[`null_two_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_two_sample.md),
[`null_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_glm.md),
[`null_mc_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_mc_fwhm.md),
[`null_mc_acf()`](https://bbuchsbaum.github.io/neurothresh/reference/null_mc_acf.md),
[`null_parcel_perm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_parcel_perm.md).

## Slots

- stat_vol:

  Observed Z-equivalent statistic map (3D array or `NeuroVol`).

- mask:

  Logical mask (3D array or `LogicalNeuroVol`).

- draw:

  `function(b)` returning a length-`sum(mask)` numeric vector for
  permutation index `1 <= b <= n_perm`. Deterministic by index.

- n_perm:

  Maximum supported draws (integer).

- df:

  Degrees of freedom (numeric, may be `NA`).

- kind:

  One of `"voxel_signflip"`, `"subject_signflip"`, `"one_sample"`,
  `"two_sample"`, `"glm"`, `"mc_fwhm"`, `"mc_acf"`, `"parcel_perm"`.

- meta:

  Free-form named list of method-specific metadata (smoothness, design
  matrix, perm indices, etc.).

# Subject-level sign-flip null generator (one-sample)

Constructs an observed Z-equivalent group map and a `null_fun(b)` that
generates sign-flip randomizations by flipping each subject map by a
single \\\pm 1\\ (constant across voxels). This preserves the spatial
covariance structure within each subject map and is the recommended way
to drive resampling-based inference for regional statistics.

## Usage

``` r
make_null_fun_subject_signflip(
  subject_vols,
  mask = NULL,
  seed = NULL,
  z_equiv = TRUE
)
```

## Arguments

- subject_vols:

  List of subject-level 3D volumes (arrays or `NeuroVol`s) with
  identical dimensions.

- mask:

  Optional logical mask aligned to the volumes. If NULL, uses all finite
  voxels common to all subjects.

- seed:

  Optional RNG seed.

- z_equiv:

  If TRUE (default), returns Z-equivalent vectors using a normal
  approximation to the one-sample t statistic with df = n_subjects - 1.

## Value

A list with components:

- z_vol:

  Observed group Z-equivalent volume (array-like, same dims).

- mask:

  Logical mask used.

- null_fun:

  Function `function(b)` returning a length-`n_mask` Z-equivalent vector
  for permutation index `b`.

- df:

  Degrees of freedom used for the t-to-Z conversion.

## Examples

``` r
subj <- replicate(10, array(rnorm(20 * 20 * 20), dim = c(20, 20, 20)), simplify = FALSE)
sf <- make_null_fun_subject_signflip(subj, seed = 1)
res <- octree_scan_fwer(sf$z_vol, mask = sf$mask, n_perm = 200, null_fun = sf$null_fun)
```

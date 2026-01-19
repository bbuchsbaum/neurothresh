# Octree scan with hierarchical step-down inference (U0-only)

Performs a multiscale dyadic-cube search using the variance-stabilized
prior-weighted mean statistic \\U_0(R)\\ and applies a Westfall-Young
step-down procedure *within each sibling family* during octree descent
(alpha-spending across the tree).

## Usage

``` r
octree_scan_stepdown(
  z_vol,
  prior_vol = NULL,
  mask = NULL,
  alpha = 0.05,
  n_perm = 1000,
  null = c("signflip_voxel", "mc_fwhm", "mc_acf"),
  null_fun = NULL,
  fwhm_vox = NULL,
  fwhm_mm = NULL,
  acf_params = NULL,
  gamma = 0.5,
  min_voxels = 8,
  min_alpha = 1e-06,
  report = c("coarsest", "finest", "greedy"),
  seed = NULL,
  two_sided = FALSE,
  prior_eta = 0.9
)
```

## Arguments

- z_vol:

  A 3D volume-like object containing Z-equivalent values.

- prior_vol:

  Optional prior weight volume aligned to `z_vol`. If NULL, uniform
  weights are used.

- mask:

  Optional logical mask (same shape as `z_vol`). If NULL, uses all
  finite voxels in `z_vol`.

- alpha:

  FWER target (default 0.05).

- n_perm:

  Number of sign-flip permutations per tested family (default 1000).

- null:

  Null calibration method.

  `"mc_fwhm"`

  :   Correlated Monte Carlo null via Gaussian smoothing of white noise
      with smoothness set by `fwhm_vox` (or `fwhm_mm` for `NeuroVol`s).
      If both are NULL, an effective `fwhm_vox` is estimated from the
      statistic map using
      [`estimate_fwhm_vox`](https://bbuchsbaum.github.io/neurothresh/reference/estimate_fwhm_vox.md).

  `"mc_acf"`

  :   Correlated Monte Carlo null using a simple mixed-ACF model
      (Gaussian + long-tailed exponential) with user-supplied parameters
      `acf_params`.

  `"signflip_voxel"`

  :   Independent voxelwise sign-flips of the group map (fast but
      generally invalid for regional/multiscale statistics; mainly
      useful for debugging).

- null_fun:

  Optional function `function(b)` returning a length `n_mask` Z-vector
  for permutation `b`. Use this to supply subject-level randomization
  that preserves spatial covariance.

- fwhm_vox:

  For `null="mc_fwhm"`, Gaussian FWHM in voxel units (length 1 or 3).

- fwhm_mm:

  For `null="mc_fwhm"`, Gaussian FWHM in mm (length 1 or 3). Requires
  `z_vol` to be a
  [`neuroim2::NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol.html)
  so voxel spacing is available.

- acf_params:

  For `null="mc_acf"`, a list with components:

  a

  :   Mixing weight in (0, 1) for the Gaussian component.

  fwhm_vox

  :   Gaussian ACF FWHM in voxel units (length 1 or 3).

  lambda_vox

  :   Exponential kernel scale in voxel units (length 1 or 3).

- gamma:

  Fraction of the local alpha budget used to test the current sibling
  family (default 0.5). Remaining budget is split among rejected
  children by prior mass.

- min_voxels:

  Minimum child size to consider for descent (default 8).

- min_alpha:

  Minimum alpha budget to continue descent (default 1e-6).

- report:

  How to select a non-redundant set of regions from the discovered
  significant nodes: `"coarsest"`, `"finest"`, or `"greedy"`.

- seed:

  Optional RNG seed.

- two_sided:

  Logical; if TRUE, uses `abs(z_vol)` for the observed map and
  permutations (two-sided via absolute Z).

- prior_eta:

  Mixing weight in \[0, 1\] to shrink the prior toward uniform mass
  (default 0.9).

## Value

A list with components:

- hits:

  Data frame of all discovered significant nodes (may include nested
  nodes).

- selected:

  Data frame of selected nodes per `report`.

- sig_mask:

  Logical mask volume for `selected` nodes.

- sig_indices:

  1-based linear indices into `z_vol` for `sig_mask`.

- params:

  Run parameters.

## Details

Compared to
[`octree_scan_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md),
this avoids a single global maxT threshold over all dyadic cubes.
Instead, it leverages the nesting of the tree and tests only small
sibling families (\<= 8) at each split.

Note: this is not the same as a global Westfall-Young step-down over the
full set of dyadic cubes (which would require the joint null
distribution over all hypotheses). Use
[`octree_scan_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
for a single-step maxT procedure over the full dyadic-cube family.

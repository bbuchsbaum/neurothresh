# Octree Max-Scan with permutation FWER (U0-only)

Computes a multiscale dyadic-cube scan statistic using the
variance-stabilized prior-weighted mean \$\$U_0(R) = \sum\_{v \in
R}\pi(v)Z(v)/\sqrt{\sum\_{v \in R}\pi(v)^2}\$\$ over all dyadic cubes in
an implicit octree, and calibrates a global threshold by a maxT
permutation (sign-flip) procedure.

## Usage

``` r
octree_scan_fwer(
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
  seed = NULL,
  two_sided = FALSE,
  prior_eta = 0.9
)
```

## Arguments

- z_vol:

  A 3D volume-like object with finite values in-brain.

- prior_vol:

  Optional prior weight volume aligned to `z_vol`. If NULL, uniform
  weights are used.

- mask:

  Optional logical mask (same shape as `z_vol`). If NULL, uses all
  finite voxels in `z_vol`.

- alpha:

  Significance level for strong FWER control (default 0.05).

- n_perm:

  Number of permutations/sign-flips (default 1000).

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

- seed:

  Optional RNG seed for reproducibility.

- two_sided:

  Logical; if TRUE, uses `abs(z_vol)` prior to sign-flips.

- prior_eta:

  Mixing weight in \[0, 1\] to shrink the prior toward uniform mass
  (default 0.9).

## Value

A list with components:

- u:

  Global maxT threshold at level `alpha`.

- M_obs:

  Observed max scan statistic.

- M_null:

  Numeric vector of permutation maxima.

- nodes:

  Data frame of all non-empty dyadic nodes with columns
  `level,i,j,k,score,scale,x0,x1,y0,y1,z0,z1`.

- sig_nodes:

  Subset of `nodes` with `score > u`.

- params:

  List of run parameters.

## Details

This implements the simplest global calibration: a single-step maxT
permutation threshold over the full dyadic-cube family.

If only a single group-level map is available, `null="mc_fwhm"` is a
safer default than voxelwise sign-flips because it preserves an assumed
spatial autocorrelation model. Exact resampling-based inference
typically requires subject-level or residual maps; supply `null_fun` to
use such randomizations.

For hierarchical step-down over sibling families during descent, see
[`octree_scan_stepdown`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_stepdown.md).

## Examples

``` r
# Map-only (model-based) calibration using a Gaussian ACF null:
z <- array(rnorm(16 * 16 * 16), dim = c(16, 16, 16))
res <- octree_scan_fwer(z, n_perm = 50, null = "mc_fwhm", fwhm_vox = c(2, 2, 2))

# User-supplied mixed ACF parameters (Gaussian + long tail):
acf <- list(a = 0.7, fwhm_vox = c(4, 4, 4), lambda_vox = c(2, 2, 2))
res2 <- octree_scan_fwer(z, n_perm = 50, null = "mc_acf", acf_params = acf)

# Subject-level sign-flips (recommended when subject maps are available):
subj <- replicate(8, array(rnorm(16 * 16 * 16), dim = c(16, 16, 16)), simplify = FALSE)
sf <- make_null_fun_subject_signflip(subj)
res3 <- octree_scan_fwer(sf$z_vol, mask = sf$mask, n_perm = 50, null_fun = sf$null_fun)
```

# Simulation Study: Resampling Choice and Null Rejection Rates

Estimates the probability of at least one false positive under a global
null for permutation-based methods, comparing:

- default voxelwise sign-flip calibration, and

- subject-level sign-flip calibration via `null_fun`.

## Usage

``` r
run_resampling_fwer_study(
  n_sims = 50,
  dims = c(12, 12, 12),
  n_subj = 20,
  n_perm = 200,
  alpha = 0.05,
  cluster_thresh = 2.5,
  fwhm_vox = 2,
  tail = c("pos", "two"),
  two_sided_policy = c("BH_all", "split_q"),
  seed = 1,
  verbose = FALSE
)
```

## Arguments

- n_sims:

  Number of independent null simulations.

- dims:

  3D volume dimensions.

- n_subj:

  Number of subject maps used to construct the group map.

- n_perm:

  Number of permutations per method call.

- alpha:

  Target level used for TFCE (`alpha`) and cluster-FDR (`q`).

- cluster_thresh:

  Cluster-forming threshold for
  [`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md).

- fwhm_vox:

  Smoothness (voxel units) used when simulating subject maps. Scalar or
  length-3 vector.

- tail:

  Tail for both methods: `"pos"` or `"two"`.

- two_sided_policy:

  Two-sided policy for
  [`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md).

- seed:

  RNG seed.

- verbose:

  If TRUE, prints simulation progress.

## Value

A list with components:

- draws:

  Per-simulation results with columns `sim`, `method`, `null_model`,
  `reject`, `time_sec`.

- summary:

  Aggregated rates and timing with columns `method`, `null_model`,
  `n_sims`, `reject_rate`, `mc_se`, `ci_low`, `ci_high`,
  `mean_time_sec`.

- params:

  Study parameters used.

## Details

Under the complete null, this "any discovery" rate coincides with the
FDR metric for
[`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md)
and is directly comparable to target level.

Subject maps are simulated as smoothed Gaussian fields and standardized
per subject. A subject-level sign-flip null is created with
[`make_null_fun_subject_signflip`](https://bbuchsbaum.github.io/neurothresh/reference/make_null_fun_subject_signflip.md)
and injected into
[`tfce_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md)
/
[`cluster_fdr`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md)
using `null_fun`.

## Examples

``` r
if (FALSE) { # \dontrun{
st <- run_resampling_fwer_study(
  n_sims = 20,
  dims = c(10, 10, 10),
  n_subj = 12,
  n_perm = 50,
  alpha = 0.1,
  seed = 1
)
st$summary
} # }
```

# Hierarchical LR-MFT Analysis

Main entry point for hierarchical likelihood-ratio matched-filter
thresholding with stepdown inference.

## Usage

``` r
hier_scan(
  stat_vol,
  prior_vol = NULL,
  parcels = NULL,
  mask = NULL,
  alpha = 0.05,
  tail = c("pos", "neg", "two"),
  kappa = c(0.5, 1, 2),
  prior_eta = 0.9,
  null_spec = NULL,
  n_perm = 1000,
  method = c("stepdown", "hierarchical"),
  gamma = 0.5,
  gamma_root = 0.5,
  min_voxels = 8,
  min_alpha = 1e-06,
  seed = NULL,
  two_sided = NULL,
  parallel = FALSE,
  n_workers = NULL
)
```

## Arguments

- stat_vol:

  A NeuroVol object containing Z-scores (or other stat type)

- prior_vol:

  Optional NeuroVol of prior weights. If NULL, uniform priors are used.

- parcels:

  Optional ClusteredNeuroVol from neuroatlas for parcel-based
  initialization. If NULL, starts from whole brain with octree descent.

- mask:

  Optional LogicalNeuroVol defining analysis domain. If NULL, uses all
  finite voxels in `stat_vol`.

- alpha:

  Significance level for FWER control (default 0.05)

- tail:

  One of `"pos"`, `"neg"`, or `"two"`.

- kappa:

  Temperature parameter(s) for soft-max. Can be a single value or vector
  of values to search over.

- prior_eta:

  Mixing weight for prior smoothing with uniform mass. Values in \[0,
  1\], where 0 is uniform and 1 is the raw prior.

- null_spec:

  Reserved; must be `NULL`. Supplied specifications are not supported
  and cause an error. Use
  [`octree_scan_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
  or
  [`octree_scan_stepdown()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_stepdown.md)
  to calibrate an octree scan with a supplied null.

- n_perm:

  Number of permutations for null distribution (default 1000)

- method:

  Inference method: "stepdown" (Westfall-Young) or "hierarchical"
  (alpha-spending only)

- gamma:

  Fraction of alpha for testing at each hierarchical level

- gamma_root:

  Fraction of alpha for root/parcel level testing (used for stepdown
  only; for hierarchical mode this is used as `gamma`)

- min_voxels:

  Minimum region size for octree descent

- min_alpha:

  Minimum alpha to continue hierarchical testing

- seed:

  Random seed for reproducibility

- two_sided:

  Deprecated logical shortcut. If TRUE and `tail = "pos"`, uses
  `tail = "two"`.

- parallel:

  Logical, whether to parallelize permutation loops using future.apply
  (optional).

- n_workers:

  Optional integer number of workers for future plan. Ignored unless
  `parallel = TRUE` and the future package is installed.

## Value

An S3 object of class "neurothresh_result" containing:

- significant_regions:

  List of significant region index vectors

- tree:

  Hierarchical tree structure with scores and p-values

- parcel_results:

  Results from parcel-level testing (if used)

- method:

  Method used

- params:

  Parameters used

## Details

The analysis proceeds as follows:

1.  Canonicalize input to Z-scores

2.  If parcels provided, test parcels with Westfall-Young step-down

3.  For significant parcels (or whole brain), descend with octree

4.  At each level, test children with step-down, spend alpha budget

5.  Collect all significant regions

## See also

[`octree_scan_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
for a global maxT scan;
[`octree_scan_stepdown`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_stepdown.md)
for octree step-down;
[`tfce_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md)
and
[`rft_peak_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/rft_peak_fwer.md)
for alternative thresholding approaches.

## Examples

``` r
if (FALSE) { # \dontrun{
library(neuroim2)
library(neuroatlas)

# Load statistical map
z_map <- read_vol("zstat1.nii.gz")

# Run with default settings
result <- hier_scan(z_map, alpha = 0.05, n_perm = 1000)

# With parcellation priors
atlas <- get_schaefer_atlas(parcels = "200", networks = "7")
result <- hier_scan(z_map, parcels = atlas$atlas)

# Summary and visualization
summary(result)
plot(result, z_map)
} # }
```

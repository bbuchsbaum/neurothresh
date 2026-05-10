# Getting Started with neurothresh

``` r

library(neurothresh)
```

## What problem are we solving?

In neuroimaging thresholding, you usually start with a statistic image:
one Z-equivalent value per voxel. A voxelwise maximum test is easy to
interpret, but it can miss spatially extended signals whose evidence is
spread across a region. `neurothresh` gives you inference tools that
keep the statistic map as the input, build a matched null distribution,
and return a standard `neurothresh_result` with a threshold, corrected
values, and a rejection mask.

The canonical workflow is:

1.  create or load a statistic map,
2.  canonicalize it to Z values if needed,
3.  choose a null specification with `null_*()`,
4.  run an inference method with `stat_vol`, `mask`, and `null_spec`,
5.  inspect the returned `neurothresh_result`.

## A first multiscale scan

Here is a small 3D example with a planted spatially extended signal. The
Monte-Carlo null is intentionally tiny so the vignette builds quickly.

``` r

set.seed(123)
z_map <- array(rnorm(20 * 20 * 20), dim = c(20, 20, 20))
z_map[8:12, 8:12, 8:12] <- z_map[8:12, 8:12, 8:12] + 2

ns <- null_mc_fwhm(
  template = z_map,
  mask = is.finite(z_map),
  smoothness = smoothness_fwhm(c(2, 2, 2)),
  n_perm = 200,
  seed = 123
)
```

``` r

fit <- octree_scan_fwer(
  stat_vol = ns$stat_vol,
  mask = ns$mask,
  null_spec = ns,
  alpha = 0.10,
  tail = "pos"
)

fit
#> <neurothresh_result>
#>   method   : octree_scan_fwer 
#>   alpha    : 0.1 
#>   tail     : pos 
#>   threshold: 4.189 
#>   rejected : 577 / 8000 voxels (7.21%)
#>   clusters : 12 
#>   n_perm   : 200
```

[`octree_scan_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
tests hierarchical spatial regions and controls familywise error over
the scanned family. The printed result tells you the method, tail,
threshold, and how many voxels are rejected at the requested alpha
level.

![The first scan starts with a statistic image, then returns a rejection
mask in the same voxel space. The white outline marks the planted signal
cube.](neurothresh_files/figure-html/octree-slice-1.png)

The first scan starts with a statistic image, then returns a rejection
mask in the same voxel space. The white outline marks the planted signal
cube.

## What is the score doing?

At each region, LR-MFT aggregates evidence across a set of voxels rather
than looking only at the single largest voxel. The lower-level scoring
helpers are useful for understanding the idea.

``` r

set.seed(42)
z_vec <- rnorm(1000)
z_vec[400:450] <- z_vec[400:450] + 2.5
pi_vec <- rep(1, length(z_vec))
region_idx <- 400:450

u0 <- score_set_stabilized(region_idx, z_vec, pi_vec)
scores <- data.frame(
  statistic = c("stabilized U0", "soft-max kappa=1", "soft-max kappa=2"),
  value = c(
    u0$U0,
    score_set(region_idx, z_vec, pi_vec, kappa = 1),
    score_set(region_idx, z_vec, pi_vec, kappa = 2)
  )
)
scores
#>          statistic     value
#> 1    stabilized U0 17.147826
#> 2 soft-max kappa=1  6.799886
#> 3 soft-max kappa=2 10.397843
```

The stabilized score accounts for region size. Increasing `kappa` makes
the soft-max score more responsive to the strongest local peaks.

## What if my map is not a Z map?

Use
[`canonicalize()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize.md)
before inference. It returns a Z-equivalent volume with the same shape
as the input.

``` r

t_map <- z_map * sqrt(25 / pmax(25 - z_map^2, 1))
z_from_t <- canonicalize(t_map, stat_type = "t", df = 25)

range(as.numeric(z_from_t), finite = TRUE)
#> [1] -4.687127  6.471432
```

[`canonicalize_stat()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize_stat.md)
is available when you need the richer object containing the Z map, mask,
and per-tail p-values. Most inference functions only need the
Z-equivalent `stat_vol`.

## Choosing a null specification

The `null_spec` object is the shared resampling seam in `neurothresh`.
It bundles an observed statistic map, a mask, and deterministic null
draws.

If you have subject-level contrast images, prefer subject-level
resampling:

``` r

ns <- null_one_sample(subject_maps, n_perm = 5000, seed = 1)

fit <- octree_scan_fwer(
  stat_vol = ns$stat_vol,
  mask = ns$mask,
  null_spec = ns
)
```

If you only have a final statistic map, use a model-based null and state
the smoothness assumption explicitly:

``` r

ns <- null_mc_fwhm(
  template = z_map,
  mask = is.finite(z_map),
  smoothness = smoothness_fwhm(8),
  n_perm = 5000,
  seed = 1
)

fit <- octree_scan_fwer(stat_vol = z_map, mask = ns$mask, null_spec = ns)
```

## How do other methods fit?

The same canonical inputs feed the comparison methods:

``` r

rft <- rft_peak_fwer(
  stat_vol = z_map,
  smoothness = smoothness_fwhm(8),
  tail = "pos"
)

tfce <- tfce_fwer(
  stat_vol = ns$stat_vol,
  mask = ns$mask,
  null_spec = ns,
  tail = "pos"
)

fdr <- cluster_fdr(
  stat_vol = z_map,
  cluster_thresh = 3.0,
  q = 0.05,
  tail = "pos",
  p_method = "rft",
  smoothness = smoothness_fwhm(8)
)
```

For distributed planned contrasts, see
[`vignette("scale-space-localized-energy")`](https://bbuchsbaum.github.io/neurothresh/articles/scale-space-localized-energy.md).
For local adaptive post-statistic filtering, see
[`vignette("lisa-local-adaptive-inference")`](https://bbuchsbaum.github.io/neurothresh/articles/lisa-local-adaptive-inference.md).

## Function map

| Task | Start with |
|----|----|
| Convert t or p maps to Z-equivalent maps | [`canonicalize()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize.md) |
| Build a subject-level one-sample null | [`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md) |
| Build a nuisance-adjusted GLM null | [`null_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_glm.md) |
| Describe field smoothness | [`smoothness_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_fwhm.md) / [`smoothness_acf()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_acf.md) |
| Run LR-MFT multiscale FWER | [`octree_scan_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md) |
| Run hierarchical stepdown scanning | [`octree_scan_stepdown()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_stepdown.md) |
| Run TFCE with a matched null | [`tfce_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md) |
| Run scale-space localized energy | [`scale_space_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md) |
| Run local adaptive LISA | [`lisa()`](https://bbuchsbaum.github.io/neurothresh/reference/lisa.md) |

# Simulation Validation: FWER and Power

## Purpose of This Vignette

Any statistical thresholding method must satisfy two requirements:

1.  **Control false positives**: When there’s no real signal (the “null
    hypothesis”), we should rarely declare significance. The Family-Wise
    Error Rate (FWER) is the probability of making even one false
    positive.

2.  **Detect real signals**: When there is a real signal, we want high
    statistical power to find it.

This vignette runs small simulations to verify that `neurothresh`
achieves both goals. We use small grids and few permutations for speed;
real analyses should use larger values.

``` r
library(neurothresh)
```

## Verifying FWER Control

First, let’s check that when there’s no signal (pure noise), we reject
at most 5% of the time when using α = 0.05. We’ll use a simple 2D
example with the soft-max score.

``` r
# Simulation parameters
B <- 200       # Number of null simulations to build threshold
n <- 30        # Grid size (30 x 30)
kernel <- make_kernel(radius = 2, sigma = 1)

# Step 1: Build null distribution by simulating pure noise
null_scores <- numeric(B)
for (b in seq_len(B)) {
  z <- matrix(rnorm(n * n), nrow = n)
  z <- smooth_field(z, kernel)  # Add spatial correlation
  null_scores[b] <- softmax_score(z, kappa = 1)
}

# Step 2: Set threshold at 95th percentile (for alpha = 0.05)
threshold <- quantile(null_scores, probs = 0.95)

# Step 3: Check false positive rate with fresh null samples
B2 <- 200
false_positives <- 0
for (b in seq_len(B2)) {
  z <- matrix(rnorm(n * n), nrow = n)
  z <- smooth_field(z, kernel)
  false_positives <- false_positives + (softmax_score(z, kappa = 1) >= threshold)
}

fwer_estimate <- false_positives / B2
cat("Estimated FWER:", fwer_estimate, "\n")
#> Estimated FWER: 0.045
cat("Expected (alpha = 0.05):", 0.05, "\n")
#> Expected (alpha = 0.05): 0.05
```

The estimated FWER should be close to 0.05. Some variation is expected
due to Monte Carlo noise.

## Power Comparison: Soft-max vs Peak Maximum

Now let’s compare two approaches for detecting a spatially extended
signal:

1.  **Soft-max score**: Aggregates evidence across the entire field
2.  **Peak maximum**: Only looks at the single highest voxel

When the signal is spread across multiple voxels (a “blob”), the
soft-max should have higher power because it accumulates evidence.

``` r
B <- 200
amp <- 1.5                    # Signal amplitude
center <- c(15, 15)           # Center of the blob

# We already have the soft-max threshold from above
thr_soft <- threshold

# Build null distribution for peak maximum
thr_max <- quantile(replicate(B, {
  z <- smooth_field(matrix(rnorm(n * n), nrow = n), kernel)
  max(z)
}), probs = 0.95)

# Now test power: how often do we detect the blob?
power_soft <- 0
power_max <- 0

for (b in seq_len(B)) {
  # Generate noise + blob signal
  z <- smooth_field(matrix(rnorm(n * n), nrow = n), kernel)
  z <- insert_blob(z, center, radius = 2, amp = amp)

  # Test both methods
  power_soft <- power_soft + (softmax_score(z, kappa = 1) >= thr_soft)
  power_max <- power_max + (max(z) >= thr_max)
}

results <- c(softmax = power_soft / B, peak_max = power_max / B)
cat("Power comparison:\n")
#> Power comparison:
cat("  Soft-max:", results["softmax"], "\n")
#>   Soft-max: NA
cat("  Peak max:", results["peak_max"], "\n")
#>   Peak max: NA
```

The soft-max typically shows higher power for this blob-shaped signal
because it pools evidence across the activated region.

## 3D Octree Scanning

The
[`octree_scan_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
function applies similar principles in 3D, using a hierarchical
decomposition. Here’s a quick check that it runs correctly:

``` r
set.seed(1)
z_null <- array(rnorm(20 * 20 * 20), dim = c(20, 20, 20))

# Run with correlated null model
res <- octree_scan_fwer(
  z_null,
  n_perm = 100,
  null = "mc_fwhm",
  fwhm_vox = c(2, 2, 2)
)

cat("Threshold:", round(res$u, 3), "\n")
#> Threshold: 7.836
cat("Observed score:", round(res$M_obs, 3), "\n")
#> Observed score: 3.81
cat("Rejected (should usually be FALSE for null data):", res$M_obs >= res$u, "\n")
#> Rejected (should usually be FALSE for null data): FALSE
```

## Practical Recommendations

These toy examples use small grids and few permutations for
demonstration. For real neuroimaging analyses:

- Use at least 1000-5000 permutations for stable thresholds
- When possible, use subject-level resampling
  (`make_null_fun_subject_signflip`) rather than model-based null
  generation
- The `mc_fwhm` null model assumes stationary Gaussian smoothness, which
  is approximate for real brain data; if you have long-tailed ACF
  parameters from upstream tools, consider `null = "mc_acf"` with
  `acf_params`

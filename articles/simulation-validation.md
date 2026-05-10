# Simulation Validation: FWER and Power

## Purpose of This Vignette

Any statistical thresholding method must satisfy two requirements:

1.  **Control false positives**: When there’s no real signal (the “null
    hypothesis”), we should rarely declare significance. The Family-Wise
    Error Rate (FWER) is the probability of making even one false
    positive.

2.  **Detect real signals**: When there is a real signal, we want high
    statistical power to find it.

This vignette runs small simulations that illustrate both goals. The
examples are deliberately small so they can run during package checks;
use larger grids, more simulation replicates, and more permutations for
a real calibration study.

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

The estimated FWER should be in the neighborhood of 0.05. Some variation
is expected because this vignette uses only 200 fresh null samples.

![The null distribution defines the decision threshold. Fresh null
simulations crossing the red line become false
positives.](simulation-validation_files/figure-html/fwer-plot-1.png)

The null distribution defines the decision threshold. Fresh null
simulations crossing the red line become false positives.

## Power Comparison: Soft-max vs Peak Maximum

Now let’s compare two approaches for detecting a spatially extended
signal:

1.  **Soft-max score**: Aggregates evidence across the entire field
2.  **Peak maximum**: Only looks at the single highest voxel

When the signal is broad and low-amplitude, the soft-max and peak
maximum can behave differently because one accumulates weak evidence and
the other responds to the single largest value.

``` r

set.seed(789)
B <- 200
amp <- 0.5                    # Low signal amplitude
center <- c(15, 15)           # Center of the blob
blob_radius <- 6              # Broad spatial extent

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
  z <- insert_blob(z, center, radius = blob_radius, amp = amp)

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

The printed values show the estimated detection rate for this particular
seeded toy setup. Change `amp` and `blob_radius` to see when aggregation
helps more than peak detection.

![Estimated detection rate for a broad low-amplitude signal. This gives
a visual check on whether region-level aggregation or peak detection is
favored under the chosen signal
shape.](simulation-validation_files/figure-html/power-plot-1.png)

Estimated detection rate for a broad low-amplitude signal. This gives a
visual check on whether region-level aggregation or peak detection is
favored under the chosen signal shape.

## 3D Octree Scanning

The
[`octree_scan_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_scan_fwer.md)
function applies similar principles in 3D, using a hierarchical
decomposition. Here’s a quick check that it runs correctly:

``` r

set.seed(1)
z_null <- array(rnorm(20 * 20 * 20), dim = c(20, 20, 20))

# Build a Monte-Carlo Gaussian null with the smoothness you want to assume
ns <- null_mc_fwhm(template = z_null, mask = is.finite(z_null),
                   smoothness = smoothness_fwhm(c(4, 4, 4)),
                   n_perm = 100, seed = 1)
res <- octree_scan_fwer(z_null, null_spec = ns)

cat("Threshold:", round(res$threshold, 3), "\n")
#> Threshold: 7.898
cat("Rejected voxels:", sum(res$reject), "\n")
#> Rejected voxels: 0
```

## Practical Recommendations

These toy examples use small grids and few permutations for
demonstration. For real neuroimaging analyses:

- Use at least 1000-5000 permutations for stable thresholds
- When possible, use subject-level resampling
  ([`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md)
  /
  [`null_subject_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_subject_signflip.md))
  rather than model-based null generation
- The `null_mc_fwhm` null model assumes stationary Gaussian smoothness,
  which is approximate for real brain data; if you have long-tailed ACF
  parameters from upstream tools, consider
  `null_mc_acf(smoothness_acf(...))` instead

## Resampling Choice for TFCE and Cluster-FDR

For TFCE and permutation cluster-FDR, the choice of resampling model can
substantially affect false-positive behavior. The helper
[`run_resampling_fwer_study()`](https://bbuchsbaum.github.io/neurothresh/reference/run_resampling_fwer_study.md)
compares:

- default voxelwise sign-flips, and
- subject-level sign-flips via
  [`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md)
  /
  [`null_subject_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_subject_signflip.md).

``` r

study <- run_resampling_fwer_study(
  n_sims = 40,
  dims = c(10, 10, 10),
  n_subj = 16,
  n_perm = 80,
  alpha = 0.05,
  cluster_thresh = 2.5,
  fwhm_vox = 2,
  tail = "pos",
  seed = 2026
)

study$summary
```

In this type of null experiment, the “any discovery” rate estimates FWER
(and equals FDR for cluster-FDR under the complete null). This gives a
practical calibration check for your resampling assumptions before
running real inference.

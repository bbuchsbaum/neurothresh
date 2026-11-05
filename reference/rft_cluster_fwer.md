# RFT Cluster-Level FWER

Computes cluster-level FWER-corrected p-values using Random Field Theory
cluster extent distribution.

## Usage

``` r
rft_cluster_fwer(
  stat_vol,
  mask = NULL,
  alpha = 0.05,
  tail = c("pos", "neg", "two"),
  cluster_thresh = 3,
  smoothness,
  df = Inf
)
```

## Arguments

- stat_vol:

  3D Z-equivalent statistic map (array or `NeuroVol`). Pass
  already-canonicalized input (use
  [`canonicalize()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize.md)
  first for t-statistics or -log10(p) maps).

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- alpha:

  Significance level (default 0.05).

- tail:

  One of `"pos"`, `"neg"`, `"two"`. For `"two"`, each direction uses
  `alpha / 2` and the union is reported, targeting FWER `alpha` across
  both tails.

- cluster_thresh:

  Cluster-forming threshold on the Z scale.

- smoothness:

  A
  [`smoothness_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_fwhm.md)
  object specifying field smoothness in millimetres. Numeric input is
  also accepted (length 1 or 3, mm) and coerced internally.

- df:

  Degrees of freedom. Use `Inf` for Gaussian field (Z-scores).

## Value

A
[`neurothresh_result`](https://bbuchsbaum.github.io/neurothresh/reference/new_neurothresh_result.md)
with `method = "rft_cluster_fwer"`. The `clusters` slot is the
per-cluster table (id, size, p_unc, p_fwer, significant). `threshold` is
`cluster_thresh`.

## Details

First thresholds at `cluster_thresh` to form clusters, then tests each
cluster's extent against the RFT null distribution. For `tail = "two"`,
the analysis runs separately on positive and negative tails at
`alpha / 2` and the union of significant clusters is reported. The
reported `p_fwer` values are doubled (capped at 1) to describe the
combined family; `p_unc` remains the one-tail cluster-extent
probability.

Extent probabilities use a three-dimensional, high-threshold Gaussian
approximation in dimensionless resels. This is an asymptotic baseline;
changing physical units does not change its probabilities.

## See also

[`rft_peak_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_peak_fwer.md),
[`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md).

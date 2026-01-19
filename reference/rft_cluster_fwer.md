# RFT Cluster-Level FWER

Computes cluster-level FWER-corrected p-values using Random Field Theory
cluster extent distribution.

## Usage

``` r
rft_cluster_fwer(
  stat_vol,
  mask = NULL,
  fwhm_mm,
  cluster_thresh = 3,
  alpha = 0.05,
  df = Inf,
  tail = c("pos", "neg", "two")
)
```

## Arguments

- stat_vol:

  NeuroVol containing statistic map

- mask:

  Optional LogicalNeuroVol for analysis domain

- fwhm_mm:

  Smoothness in mm

- cluster_thresh:

  Cluster-forming threshold (Z-score)

- alpha:

  Significance level (default 0.05)

- df:

  Degrees of freedom. Use Inf for Gaussian field.

- tail:

  "pos", "neg", or "two"

## Value

List with components:

- sig_mask:

  Logical mask of significant clusters

- table:

  Data frame with cluster statistics

- n_clusters:

  Number of significant clusters

## Details

First thresholds at cluster_thresh to form clusters, then tests each
cluster's extent against the RFT null distribution.

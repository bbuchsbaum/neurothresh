# Cluster-FDR with RFT or Permutation P-values

Controls the False Discovery Rate at the cluster level by applying the
Benjamini-Hochberg procedure to cluster-extent p-values.

## Usage

``` r
cluster_fdr(
  stat_vol,
  mask = NULL,
  cluster_thresh = 3,
  q = 0.05,
  p_method = c("rft", "perm"),
  fwhm_mm = NULL,
  n_perm = 5000,
  tail = c("pos", "neg", "two"),
  two_sided_policy = c("BH_all", "split_q"),
  seed = NULL
)
```

## Arguments

- stat_vol:

  NeuroVol containing statistic map

- mask:

  Optional LogicalNeuroVol for analysis domain

- cluster_thresh:

  Cluster-forming threshold (Z-score)

- q:

  FDR level (default 0.05)

- p_method:

  Method for computing cluster p-values: "rft" (Random Field Theory) or
  "perm" (permutation)

- fwhm_mm:

  Smoothness in mm (required for p_method = "rft")

- n_perm:

  Number of permutations (for p_method = "perm")

- tail:

  "pos", "neg", or "two"

- two_sided_policy:

  For two-sided: "BH_all" (pool both tails) or "split_q" (q/2 per tail)

- seed:

  Random seed for permutation

## Value

List with components:

- sig_mask:

  Logical mask of FDR-significant clusters

- table:

  Data frame with cluster statistics

- n_clusters:

  Number of significant clusters

## Details

This implements "topological FDR" as advocated by Chumbley & Friston,
where FDR is controlled over the finite set of clusters rather than over
all voxels.

## References

Chumbley, J. R., & Friston, K. J. (2009). False discovery rate
revisited: FDR and topological inference using Gaussian random fields.
NeuroImage, 44(1), 62-70.

## Examples

``` r
if (FALSE) { # \dontrun{
# RFT-based cluster FDR
result <- cluster_fdr(z_map, cluster_thresh = 3.0, q = 0.05,
                      p_method = "rft", fwhm_mm = 8)

# Permutation-based cluster FDR
result <- cluster_fdr(z_map, cluster_thresh = 3.0, q = 0.05,
                      p_method = "perm", n_perm = 5000)
} # }
```

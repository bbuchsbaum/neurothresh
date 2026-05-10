# Cluster-FDR with RFT or Permutation P-values

Controls the False Discovery Rate at the cluster level by applying the
Benjamini-Hochberg procedure to cluster-extent p-values.

## Usage

``` r
cluster_fdr(
  stat_vol,
  mask = NULL,
  q = 0.05,
  tail = c("pos", "neg", "two"),
  cluster_thresh = 3,
  p_method = c("rft", "perm"),
  smoothness = NULL,
  null_spec = NULL,
  n_perm = NULL,
  seed = NULL,
  fwhm_mm = NULL,
  null_fun = NULL,
  allow_voxel_signflip = FALSE,
  two_sided_policy = c("BH_all", "split_q")
)
```

## Arguments

- stat_vol:

  3D Z-equivalent statistic map (array or `NeuroVol`).

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- q:

  FDR level (default 0.05).

- tail:

  One of `"pos"`, `"neg"`, `"two"`. For `"two"`, both tails are pooled
  into a single BH ordering.

- cluster_thresh:

  Cluster-forming threshold (Z-score).

- p_method:

  Method for computing cluster p-values: `"rft"` (Random Field Theory)
  or `"perm"` (permutation).

- smoothness:

  A
  [`smoothness_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_fwhm.md)
  (or numeric mm) — required when `p_method = "rft"`. Ignored otherwise.

- null_spec:

  A
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  supplying null draws — required when `p_method = "perm"`.

- n_perm:

  Number of permutations to consume from `null_spec`. Defaults to
  `null_spec$n_perm`.

- seed:

  Optional random seed.

- fwhm_mm:

  Deprecated alias for `smoothness`.

- null_fun:

  Deprecated bare-function null source; prefer `null_spec`.

- allow_voxel_signflip:

  Deprecated; if `TRUE` and `null_spec`/`null_fun` are missing, falls
  back to independent voxel sign-flips.

- two_sided_policy:

  Deprecated. The `"two"` tail now uses BH on the pooled set of clusters
  (the previous `"BH_all"` behaviour). Pass `"split_q"` only if you need
  the legacy q/2 split per tail.

## Value

A
[`neurothresh_result`](https://bbuchsbaum.github.io/neurothresh/reference/new_neurothresh_result.md)
with `method = "cluster_fdr"` and the cluster table in `$clusters`.
`$threshold` is `cluster_thresh`.

## References

Chumbley & Friston (2009), NeuroImage 44(1), 62-70.

## See also

[`rft_cluster_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_cluster_fwer.md),
[`tfce_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md),
[`hier_scan()`](https://bbuchsbaum.github.io/neurothresh/reference/hier_scan.md).

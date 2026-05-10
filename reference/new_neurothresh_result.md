# Construct a neurothresh result

Internal constructor used by every public inference method to return a
uniform result object. The S3 class is `"neurothresh_result"`;
downstream code should treat the object as opaque except for the
documented slots.

## Usage

``` r
new_neurothresh_result(
  method,
  stat_vol,
  mask,
  threshold,
  reject,
  pvals,
  clusters = NULL,
  null_max = NULL,
  params = list(),
  call = NULL,
  extras = list()
)
```

## Arguments

- method:

  Character; one of `"hier_scan"`, `"octree_scan_fwer"`,
  `"octree_scan_stepdown"`, `"tfce_fwer"`, `"cluster_fdr"`,
  `"rft_peak_fwer"`, `"rft_cluster_fwer"`, or `"lisa"`.

- stat_vol:

  The observed (canonicalized) statistic map.

- mask:

  Logical mask used.

- threshold:

  Critical value on the relevant statistic (numeric scalar).

- reject:

  Logical 3D rejection volume.

- pvals:

  Numeric 3D corrected p-value volume.

- clusters:

  Optional data frame of clusters (NULL for voxelwise-only methods).

- null_max:

  Optional vector of max-null statistics for permutation methods.

- params:

  List of method-specific parameters (alpha/q, tail, n_perm, kappa, H/E,
  smoothness, etc.).

- call:

  The captured [`match.call()`](https://rdrr.io/r/base/match.call.html)
  from the public method.

- extras:

  Optional named list of method-specific extras (e.g., per-node tables
  for hier_scan, TFCE map for tfce_fwer). Stored under `$extras`.

## Value

Object of class `"neurothresh_result"`.

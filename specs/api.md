# Public API Specification

## Main Entry Points

### hier_scan()
Primary function for hierarchical LR-MFT analysis.

```r
hier_scan(
  z_vol,              # NeuroVol: Z-score map
  prior_vol = NULL,   # NeuroVol: prior weights (NULL = uniform)
  parcels = NULL,     # ClusteredNeuroVol: parcellation (NULL = octree only)
  mask = NULL,        # LogicalNeuroVol: analysis mask (NULL = non-zero voxels)
  alpha = 0.05,       # Numeric: FWER level
  kappa = 1.0,        # Numeric: temperature parameter
  n_perm = 1000,      # Integer: permutations (0 = parametric)
  method = "stepdown",# Character: "stepdown" or "hierarchical"
  min_voxels = 8      # Integer: minimum region size for octree
)
```

**Returns:** S3 object of class `"neurothresh_result"`

---

## Scoring Functions

### score_set()
```r
score_set(indices, z_vol, prior_vol, kappa = 1.0)
```
Returns: T_κ(R) = log Σ π·exp(κ·Z)

### score_set_stabilized()
```r
score_set_stabilized(indices, z_vol, prior_vol)
```
Returns: U_0(R) = Σπ·Z / sqrt(Σπ²)

---

## Inference Functions

### generate_null_scores()
```r
generate_null_scores(z_vol, prior_vol, regions, n_perm, kappa = 1.0, seed = NULL)
```
Returns: Matrix [n_perm × n_regions]

### wy_stepdown()
```r
wy_stepdown(observed, null_matrix, alpha = 0.05)
```
Returns: data.frame(region, score, p_adj, rejected)

### hier_test()
```r
hier_test(node, alpha, null_dist, method = "bonferroni")
```
Returns: Tree with rejection status

---

## Baseline Methods

### RFT
```r
rft_peak_fwer(z_vol, mask, fwhm, alpha = 0.05)
rft_cluster_fwer(z_vol, mask, fwhm, cluster_thresh, alpha = 0.05)
```

### TFCE
```r
tfce_transform(z_vol, mask, H = 2.0, E = 0.5, dh = 0.1)
tfce_fwer(z_vol, mask, n_perm = 5000, alpha = 0.05)
```

### Cluster-FDR
```r
cluster_fdr(z_vol, mask, cluster_thresh, q = 0.05, p_method = "rft")
```

---

## Utility Functions

### canonicalize_stat()
```r
canonicalize_stat(vol, type = c("Z", "t", "neglog10p"), df = NULL)
```
Converts any stat type to Z-scores.

### load_atlas_priors()
```r
load_atlas_priors(atlas)
```
Extracts ClusteredNeuroVol from neuroatlas object.

---

## Two-Sided Variants

All inference functions have `_twosided` variants:
```r
hier_scan_twosided(...)
rft_peak_fwer_twosided(...)
tfce_fwer_twosided(...)
```

---

## Result Object

Class `"neurothresh_result"` contains:
- `$significant_regions` - List of significant region index vectors
- `$tree` - Hierarchical tree structure with scores/p-values
- `$threshold` - Effective threshold used
- `$method` - Method used ("lrmft", "rft", "tfce", "cluster_fdr")
- `$params` - Parameters used (alpha, kappa, n_perm, etc.)

Methods:
- `print.neurothresh_result(x)`
- `summary.neurothresh_result(x)`
- `plot.neurothresh_result(x, vol, ...)`
- `as_mask.neurothresh_result(x, space)` - Returns LogicalNeuroVol

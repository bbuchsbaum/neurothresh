# neurothresh NEWS

## neurothresh 0.2.0 (unreleased) — API regularization

This release regularises every public inference method behind a single
canonical signature. The plan and motivation live in
`specs/api-consistency-plan.md`.

### One stat-map name

Every inference method takes its observed map as the **first positional
argument**, named `stat_vol`. The old name `z_vol` is gone.

```r
# Before
hier_scan(z_vol = z_map, ...)
octree_scan_fwer(z_vol = z_map, ...)

# After
hier_scan(stat_vol = z_map, ...)        # or just hier_scan(z_map, ...)
octree_scan_fwer(stat_vol = z_map, ...)
```

### One tail vocabulary

`tail = c("pos", "neg", "two")` everywhere. The `two_sided` boolean is
deprecated; passing `two_sided = TRUE` still works but is mapped onto
`tail = "two"`.

`rft_fwer_twosided()` is removed — pass `tail = "two"` to `rft_peak_fwer()`
or `rft_cluster_fwer()` instead. The `two_sided_policy` argument on
`cluster_fdr()` is also deprecated; the default `"BH_all"` semantics now
follows from `tail = "two"`.

### One null seam

A new S3 class `neurothresh_null_spec` is the uniform handle every
permutation-based method consumes. Build it with one of the design-named
factories:

```r
ns <- null_one_sample(subject_maps, n_perm = 5000)
ns <- null_two_sample(subject_maps, group, n_perm = 5000)
ns <- null_glm(subject_maps, design = ~ group + age, contrast = "groupcase")
ns <- null_voxel_signflip(z_map, n_perm = 5000)
ns <- null_mc_fwhm(z_map, mask, smoothness = smoothness_fwhm(8), n_perm = 5000)
ns <- null_mc_acf(z_map, mask, smoothness = smoothness_acf(0.7, 3, 8), n_perm = 5000)
ns <- null_parcel_perm(z_map, parcels, n_perm = 5000)
ns <- null_subject_signflip(subject_maps, n_perm = 5000)   # alias of null_one_sample
```

Inference methods now accept `null_spec` as their canonical null source.
The deterministic-by-index draw function lives at `ns$draw`.

The `make_null_fun_*()` factories are kept as deprecated aliases for one
release; they internally delegate to the new factories. The lazy-signs
variant of `make_null_fun_subject_signflip()` has been retired — its
replacement uses the deterministic-by-index sign matrix from
`null_one_sample()`.

```r
# Recommended migration for TFCE
ns <- null_one_sample(subject_maps, n_perm = 5000, seed = 1)
tfce <- tfce_fwer(ns$stat_vol, mask = ns$mask, null_spec = ns, tail = "pos")

# Recommended migration for cluster-FDR (perm)
ns <- null_one_sample(subject_maps, n_perm = 5000, seed = 1)
fdr <- cluster_fdr(ns$stat_vol, mask = ns$mask, cluster_thresh = 3.0,
                   q = 0.05, p_method = "perm", null_spec = ns, tail = "pos")
```

### One smoothness type

A new S3 class `neurothresh_smoothness` carries field-smoothness
specifications:

```r
smoothness_fwhm(8)              # length-1 isotropic mm
smoothness_fwhm(c(6, 6, 8))     # per-axis mm
smoothness_acf(a = 0.7, b = 3.0, c = 8.0)
```

Methods that need smoothness now accept these objects through a single
`smoothness =` argument. The old positional / keyword variants
(`fwhm_mm =`, `fwhm_vox =`, `acf_params =`) remain for one release as
deprecated aliases.

```r
# Before
rft_peak_fwer(z_map, fwhm_mm = 8)
cluster_fdr(z_map, p_method = "rft", fwhm_mm = 8, ...)

# After
rft_peak_fwer(z_map, smoothness = smoothness_fwhm(8))
cluster_fdr(z_map, p_method = "rft", smoothness = smoothness_fwhm(8), ...)
```

### Canonicalize is the user's first step

`canonicalize()` is now a top-level exported function returning the
Z-equivalent volume directly:

```r
zvol <- canonicalize(t_map, stat_type = "t", df = 24)
res <- hier_scan(zvol, parcels = atlas, alpha = 0.05)
```

`canonicalize_stat()` is still exported for callers who want the richer
list (one-/two-sided p-values plus mask). `hier_scan()` no longer accepts
`stat_type` or `df` — pass an already-canonicalized map.

### One result class

Every public method returns an S3 `neurothresh_result` with consistent
slots: `method`, `stat_vol`, `mask`, `threshold`, `reject`, `pvals`,
`clusters`, `null_max`, `params`, `call`, `extras`. `print()`,
`summary()`, and `as.data.frame()` methods are provided. Old return
shapes (`$sig_mask`, `$u`, `$M_obs`, `$max_null`, `$table`,
`$n_clusters`, `$cluster_labels`, `$tfce`) remain available as aliases on
the new object so existing code keeps working.

### Scale-space frequency diagnostics

Scale-space inference now separates localization evidence from spatial-frequency
readout. `local_energy()` remains the localized squared-evidence channel.
New `global_lowpass_energy()` and `global_band_energy()` channels support
cumulative low-pass and band-pass frequency curves, respectively.
`scale_space_diagnostics()` labels each requested curve so users can tell
whether a result should be interpreted as localization evidence or as a
spatial-frequency readout.

### LISA uses one bilateral filter path

LISA now uses a single production filter implementation:
`neuroim2::bilateral_filter()` from neuroim2 0.16.0, whose bilateral
neighborhoods are mask-normalized by default. The public `backend` switch and
reserved sparse backend have been removed. The filter passes a fixed
`range_scale` when using the default `scale = "sd_first_30_nulls"` policy, so
observed and null statistic maps are filtered under the same intensity
bandwidth.

### Polymorphic subject-level inputs

`as_subject_matrix()` (replacing `as_image_matrix()`) now accepts a list
of `NeuroVol`s, a 4D array, a `NeuroVec`, a `SparseNeuroVec`, or a
subject-by-voxel matrix. All subject-level null factories
(`null_one_sample()`, `null_two_sample()`, `null_glm()`,
`null_subject_signflip()`) route through it.

Mask precedence is strict: when a mask is supplied, any non-finite
voxel inside it is an error, and for `SparseNeuroVec` any voxel outside
the intrinsic mask is an error.

### Porting recipe in one block

```r
# Old
make_null_fun_subject_signflip(subj)                        # → null_one_sample(subj, n_perm = 5000)
make_null_fun_voxel_signflip(z[mask])                       # → null_voxel_signflip(z, mask)
make_null_fun_mc_fwhm(..., fwhm_vox = fv)                   # → null_mc_fwhm(template, mask, smoothness_fwhm(...))
make_null_fun_mc_acf(..., acf_params = ap)                  # → null_mc_acf(template, mask, smoothness_acf(...))
canonicalize_stat(z, stat_type = "t", df = 24)$Zeq          # → canonicalize(z, stat_type = "t", df = 24)
rft_peak_fwer(z, fwhm_mm = 8)                               # → rft_peak_fwer(z, smoothness = smoothness_fwhm(8))
rft_fwer_twosided(z, method = "peak", fwhm_mm = 8)          # → rft_peak_fwer(z, tail = "two", smoothness = smoothness_fwhm(8))
tfce_fwer(z, allow_voxel_signflip = TRUE)                   # → tfce_fwer(z, null_spec = null_voxel_signflip(z))
hier_scan(z, two_sided = TRUE)                              # → hier_scan(z, tail = "two")
hier_scan(z, stat_type = "t", df = 24)                      # → hier_scan(canonicalize(z, "t", 24))
```


## neurothresh 0.1.1

### Breaking change (February 22, 2026)

- `tfce_fwer()` now requires `null_fun` by default.
- `cluster_fdr(..., p_method = "perm")` now requires `null_fun` by default.
- Unsafe voxelwise sign-flip calibration is still available, but only with explicit opt-in:
  - `allow_voxel_signflip = TRUE`

This change hardens inference defaults and avoids accidentally running anti-conservative permutation calibration for regional statistics.

### Migration guide

Recommended migration for TFCE:

```r
null_gen <- make_null_fun_subject_signflip(subject_maps, seed = 1)

res_tfce <- tfce_fwer(
  null_gen$z_vol,
  mask = null_gen$mask,
  n_perm = 5000,
  null_fun = null_gen$null_fun,
  tail = "pos"
)
```

Recommended migration for permutation cluster-FDR:

```r
null_gen <- make_null_fun_subject_signflip(subject_maps, seed = 1)

res_fdr <- cluster_fdr(
  null_gen$z_vol,
  mask = null_gen$mask,
  cluster_thresh = 3.0,
  q = 0.05,
  p_method = "perm",
  n_perm = 5000,
  null_fun = null_gen$null_fun,
  tail = "pos"
)
```

If you must preserve previous behavior temporarily:

```r
tfce_fwer(z_map, n_perm = 5000, allow_voxel_signflip = TRUE)

cluster_fdr(
  z_map,
  cluster_thresh = 3.0,
  q = 0.05,
  p_method = "perm",
  n_perm = 5000,
  allow_voxel_signflip = TRUE
)
```

# Scale-Space Energy Inference Notes

## Assessment

This method is a strong fit for `neurothresh`. The package already compares regional LR-MFT, octree max scans, TFCE, RFT, cluster-FDR, and permutation-calibrated FWER. Scale-space localized energy inference fills a real gap: atlas-free, resolution-aware evidence for distributed effects. It should be presented as a separate method family, not as a replacement for `hier_scan()`: LR-MFT asks which hierarchical regions carry evidence, while scale-space localized energy asks the finest spatial resolution at which an effect is still detectable after correction over scale and location.

The best package framing is:

> Scale-space localized energy inference for distributed interaction maps.

The central output should be a resolution-evidence summary:

- omnibus scale-space `p_FWER`
- best-evidence scale
- finest detected scale
- stable supported scale range
- scale-specific corrected maps
- energy-mass summaries at the supported scale

## Statistical Target

The input can be either a precomputed statistic map plus calibrated null generator, or subject-level contrast images. For the group x task x model interaction, start from contrast images `D_i(v)`, then fit a group model:

```text
D_i(v) = alpha(v) + beta_G(v) G_i + C_i gamma(v) + epsilon_i(v)
```

The method scans the resulting statistic map `z(v)` across heat-kernel or Gaussian smoothing scales. Large scales test coarse distributed effects; small scales approach voxelwise localization. This first version should avoid true spatial-frequency claims unless a later band-pass or graph-wavelet channel is added.

## Proposed User API

Low-level, package-consistent entry point:

```r
fit <- scale_space_energy(
  stat_vol,
  mask = gray_matter_mask,
  scales = scale_space(
    fwhm_mm = c(6, 10, 14, 18, 24, 32, 48),
    include_identity = TRUE,
    include_global = TRUE
  ),
  statistics = list(
    local_energy(),
    signed_coherence(tail = "two")
  ),
  null_fun = null_fun,
  n_perm = 5000,
  correction = "max_over_scale_location_statistic"
)
```

High-level follow-on:

```r
fit <- scale_space_glm(
  contrast_images = D_task_model,
  group = group,
  covariates = ~ age + sex + motion,
  permutation = permutation_glm(method = "freedman_lane", n_perm = 10000)
)
```

Start with the low-level API because it matches current `tfce_fwer()` and `octree_scan_fwer()` contracts.

## Core Algorithm

1. Build a graph over mask voxels using physical voxel coordinates.
2. Form a graph Laplacian `L`.
3. For each scale, apply `S_rho = exp(-t(rho) L)` to the statistic vector.
4. Compute channels:
   - `signed_coherence`: `S_rho z`
   - `local_energy`: `S_rho z^2`
   - `positive_energy`: `S_rho pmax(q * z, 0)^2`
5. Repeat for every permutation draw from `null_fun`.
6. Store permutation maxima over scale, location, and channel.
7. Convert observed channel maps to corrected p-value maps using the max-null distribution.
8. Report best scale, finest detected scale, and stable supported scale range.

Important parameterization detail: if user-facing `rho` is FWHM in mm, convert to Gaussian sigma with `sigma = fwhm / sqrt(8 * log(2))`. The graph heat time `t` must be calibrated to the chosen Laplacian and voxel spacing; do not silently assume `t = rho^2 / 2` unless tests show the graph impulse response matches the requested FWHM. Include the global energy endpoint and identity endpoint from the start.

## Test Statistics

Raw energy fields are easy to explain but not always comparable across scales and locations. The implementation should support studentized channel values:

```text
T_energy(rho, v) = [E_rho(v) - mu_null(rho, v)] / sd_null(rho, v)
```

The null mean and variance can be estimated from permutations. This is especially important when correcting over mixed channels. Raw maps can still be returned for interpretation and energy-mass plots.

For local energy, report raw, excess, and studentized maps. Energy-mass summaries should use above-null energy, such as `excess_energy` or the positive part of the studentized map, not raw `S_rho z^2`.

## Reuse From `adjoin`

Useful pieces to borrow or depend on:

- `spatial_adjacency()` for sparse coordinate-based adjacency.
- `spatial_laplacian()` and `laplacian.neighbor_graph()` for graph Laplacians.
- `graph_weights_fast()` for scalable kNN graph construction.
- `normalize_adjacency()` and isolate handling patterns.

Do not reuse `compute_diffusion_kernel()` directly for whole-brain inference: it returns a dense kernel matrix and is not suitable for tens of thousands of voxels. Prefer one of:

- spectral truncation with `RSpectra` and cached eigenpairs;
- Krylov/Chebyshev heat-kernel vector products;
- separable masked Gaussian smoothing as a fast regular-grid backend.

`adjoin` is best treated as a source of graph-building utilities, not the final smoothing engine.

## Integration With Existing Package

The method should reuse current helper concepts:

- `null_fun(b)` returning a mask-space statistic vector;
- `.coerce_null_draw()` style conversion to mask-space vectors;
- `fwer_threshold_from_maxnull()` and `fwer_p_from_maxnull()` for discrete permutation calibration;
- `.array_to_vol()` for returning volume-like maps;
- subject-level null helpers for one-sample contrasts.

Needed new support:

- `make_null_fun_glm()` or `permutation_glm()` for nuisance-adjusted designs;
- `make_within_subject_contrast_images()` for repeated-measures factorial contrast images;
- a mask-to-coordinate graph builder for `neuroim2`/array inputs;
- a scale-space object recording FWHM, sigma, effective basis count, and backend;
- plotting helpers for resolution curves and supported-scale maps.

## Outputs

`summary.scale_space_energy()` should report:

- omnibus corrected p-value;
- strongest channel and scale;
- finest detected scale;
- stable supported scale range;
- effective spatial basis count at each scale;
- whether voxelwise/TFCE comparators were significant if supplied.

Map outputs should distinguish inferential maps from descriptive energy-mass maps. A corrected point `(rho, v, channel)` means a localized test centered at `v` with bandwidth `rho` is significant; it does not mean the individual voxel is significant.

## Validation Plan

Minimum tests before exposing the API:

- `rho = 0` or identity backend reproduces voxelwise signed and energy endpoints.
- Increasing FWHM decreases effective basis count.
- Tiny graphs match dense `expm` results within tolerance.
- Spectral and approximate backends agree on small masks.
- Corrected p-values follow the plus-one permutation rule.
- Null simulations are not anti-conservative for one-sample sign-flip maps.
- Local energy detects mixed-sign distributed blobs better than signed coherence.
- Signed coherence detects coherent same-sign blobs.

## Recommended Implementation Order

1. Implement a small internal Gaussian scale-space smoother on arrays.
2. Add `scale_space_energy(stat_vol, mask, null_fun, ...)`.
3. Include identity and global endpoints.
4. Add two-pass studentization, summary, and resolution-curve plotting.
5. Add energy-mass descriptive maps at selected scales.
6. Add one-sample/two-group wrappers and repeated-measures contrast-image helpers.
7. Add `make_null_fun_glm()` for group/covariate models.
8. Add bootstrap stability as a downstream descriptive layer.

See `specs/scale-space-energy-plan.md` for the full modular implementation
plan, including one-sample and two-group raw-image resampling support.

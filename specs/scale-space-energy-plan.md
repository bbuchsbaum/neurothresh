# Scale-Space Localized Energy Inference Implementation Plan

## Goal

Add an elegant, modular method family for atlas-free distributed-effect inference:

> Scale-space localized energy inference for distributed interaction maps.

The method should fit `neurothresh` as a permutation-calibrated bridge between global energy tests and voxelwise localization. It should answer: at what spatial resolution is this planned effect detectable with familywise-corrected evidence?

This is not a replacement for LR-MFT or octree inference. It is a companion method for cases where effects are expected to be spatially distributed, possibly mixed-sign, and not well represented by an atlas or a hard cluster threshold.

Avoid claiming that the first version estimates an exact spatial frequency. Gaussian scale-space smoothing is a low-pass resolution/localization family: effects visible at one FWHM can also contribute to coarser smoothed maps. The first-version inferential claim is therefore "detectable at this spatial resolution scale." A later band-pass or graph-wavelet channel can support stronger spatial-frequency claims.

## Design Principles

Keep the implementation layered:

1. Observation/model layer: raw participant images to effect statistic maps.
2. Resampling layer: synchronized null statistic maps.
3. Scale-space layer: smoothing operators over graph/grid scales.
4. Statistic-channel layer: signed coherence, local energy, positive energy.
5. Correction/output layer: max-null FWER over scale, location, and channel.

The low-level scan should only require `stat_vol`, `mask`, and `null_fun`, matching the package's current `tfce_fwer()` and `octree_scan_fwer()` style. Higher-level raw-image design helpers should produce those inputs.

This mirrors the useful part of the LISA precedent: apply spatial filtering after the voxelwise statistic map, apply the same operation to permutation maps, and infer on the transformed maps. The difference is that this method searches a prespecified scale-space family and controls FWER over scale x location x channel, whereas LISA uses a post-statistic adaptive filter and controls FDR in filtered maps.

## Public API Shape

### Low-Level Scanner

```r
fit <- scale_space_energy(
  stat_vol,
  mask = NULL,
  scales = scale_space(
    fwhm_mm = c(6, 10, 14, 18, 24, 32, 48),
    include_identity = TRUE,
    include_global = TRUE
  ),
  statistics = list(
    local_energy(),
    signed_coherence(tail = "two")
  ),
  null_fun,
  n_perm = 5000,
  studentize = TRUE,
  backend = c("grid_gaussian", "graph_heat", "spectral"),
  alpha = 0.05,
  seed = NULL
)
```

This is the first implementation target. It is design-agnostic and can be used with one-sample, two-group, or externally supplied GLM null maps.

### One-Sample Raw-Image Helper

```r
ss <- scale_space_onesample(
  subject_images,
  mask = NULL,
  scales = scale_space(fwhm_mm = exp(seq(log(4), log(64), length.out = 24))),
  statistics = list(local_energy(), signed_coherence(tail = "two")),
  n_perm = 5000,
  resampling = one_sample_signflip()
)
```

This should wrap a stronger version of `make_null_fun_subject_signflip()`: build an observed one-sample t/Z map from subject images, pre-generate subject-level sign flips, and return reproducible null maps.

### Two-Group Raw-Image Helper

```r
ss <- scale_space_twosample(
  subject_images,
  group,
  mask = NULL,
  covariates = NULL,
  scales = scale_space(fwhm_mm = c(6, 10, 14, 18, 24, 32, 48)),
  n_perm = 5000,
  resampling = two_sample_label_permutation()
)
```

For the first two-group version, support unpaired two-level group labels without nuisance covariates. With covariates, unbalanced clinical/control data, or plausible group variance differences, route through the GLM permutation helper below. Simple label permutation is a convenience wrapper, not the preferred SDAM/control analysis path.

### General GLM Helper

```r
glm_null <- make_null_fun_glm(
  images = subject_images,
  design = ~ group + age + sex + motion,
  contrast = "groupcase",
  mask = gray_matter_mask,
  permutation = permutation_glm(method = "freedman_lane", n_perm = 5000),
  statistic = c("t", "z")
)

fit <- scale_space_energy(glm_null$stat_vol, mask = glm_null$mask,
                          null_fun = glm_null$null_fun, n_perm = 5000)
```

This is the right extension point for group differences with nuisance covariates and planned interaction contrasts.

### Repeated-Measures Contrast Helper

```r
D_task_model <- make_within_subject_contrast_images(
  images,
  subject = "subject_id",
  within = c("task", "model"),
  contrast = "task:model"
)
```

This helper should turn repeated-measures factorial image data into one planned contrast image per subject. The scale-space machinery then sees only the resulting planned interaction contrast map or the group-level null generator built from those contrast images.

## Data Contracts

Create one internal representation for image observations:

```r
as_image_matrix(images, mask) -> list(
  X = numeric matrix [n_subjects x n_mask],
  mask = logical volume,
  mask_idx = integer,
  template = original first image,
  dims = integer length 3,
  coords = numeric matrix [n_mask x 3]
)
```

Accept lists of arrays, lists of `NeuroVol`s, or an already prepared numeric matrix plus a template/mask. The raw-image helpers should use this contract, while `scale_space_energy()` should keep accepting a precomputed map plus `null_fun`.

## Resampling Design

### One-Sample

Use subject-level sign flips. Each permutation multiplies a subject's entire image by `+1` or `-1`; voxelwise independent sign flips should remain only an explicit unsafe fallback in low-level testing code.

Observed statistic:

```text
t(v) = mean(D_i(v)) / [sd(D_i(v)) / sqrt(n)]
```

Then convert to Z-equivalent if requested. The null function recomputes the same statistic after sign flips.

### Two-Group

For simple unpaired two-group designs, use label permutation with the group sizes fixed. The statistic should be a two-sample t map, with Welch and pooled-variance options:

```r
two_sample_label_permutation(variance = c("welch", "pooled"))
```

Default to Welch unless pooled variance is explicitly requested. Pre-generate permutation indices so `null_fun(b)` is deterministic.

Documentation should state that label permutation assumes exchangeability under the null. For unbalanced clinical comparisons or nuisance-adjusted models, prefer `make_null_fun_glm(..., permutation = permutation_glm(method = "freedman_lane"))`. Add a roadmap item for wild/bootstrap-style or heteroskedasticity-robust residual resampling.

### GLM With Nuisance Covariates

Implement Freedman-Lane first:

1. Fit the reduced model without the tested effect.
2. Permute residuals under exchangeability constraints.
3. Add permuted residuals to reduced fitted values.
4. Refit the full model and compute the contrast statistic map.

Expose the resampling object, but keep the low-level scanner independent of it. Later methods can include ter Braak or sign-flip residual methods for symmetric errors.

### Exchangeability Blocks

Reserve API space for:

```r
permutation_glm(..., block = subject_family, strata = site)
```

Do not implement complex blocking in the first pass, but avoid APIs that make it impossible.

## Scale-Space Operators

Provide a `scale_space()` constructor:

```r
scale_space(
  fwhm_mm,
  include_identity = TRUE,
  include_global = TRUE,
  report = c("fwhm", "effective_basis")
)
```

Each scale records FWHM, Gaussian sigma, backend time parameter, and effective basis count.

The global endpoint is part of the initial conceptual contract, not an optional later embellishment. It makes the continuum explicit:

```text
global energy -> coarse local energy -> medium local energy -> fine local energy -> identity/voxelwise endpoint
```

### Backend Order

1. `grid_gaussian`: first backend for regular 3D volumes. Fast, simple, and aligned with user-facing FWHM.
2. `graph_heat`: graph Laplacian backend over the mask for irregular domains.
3. `spectral`: cached truncated eigenbasis for effective-basis reporting and repeated scans.

Borrow from `~/code/adjoin` for sparse adjacency, Laplacian construction, kNN graph utilities, and isolate handling. Do not depend on dense diffusion kernels for full images.

## Statistic Channels

Implement channels as independent objects/functions with a common signature:

```r
channel_fun(z_vec, smoother, direction) -> numeric vector length n_mask
```

Initial channels:

- `signed_coherence`: `S_rho z`
- `local_energy`: `S_rho z^2`
- `positive_energy`: `S_rho pmax(q * z, 0)^2`

Use `local_energy` as the primary channel for mixed-sign distributed patterns. Preserve signed coherence because same-sign distributed effects should not be forced into squared-energy interpretation.

Tail handling must be channel-specific. A global `tail = "neg"` is meaningful for signed coherence but ambiguous for local energy, which is already unsigned. Use channel constructors such as:

```r
statistics = list(
  signed_coherence(tail = "two"),
  local_energy(),
  positive_energy(direction = "control_gt_sdam")
)
```

`positive_energy()` should require explicit prespecification because it adds another searched channel.

Reserve a later true spatial-frequency channel using a difference-of-Gaussians or graph-wavelet band-pass:

```text
B_rho z = S_{rho / sqrt(2)} z - S_{rho * sqrt(2)} z
E_band_rho(v) = S_{eta rho}[(B_rho z)^2](v)
```

## Studentization

Default `studentize = TRUE` for inferential statistics:

```text
T(rho, v, c) = [A_obs(rho, v, c) - mu_null(rho, v, c)] / sd_null(rho, v, c)
```

This makes channels and scales more comparable before taking the max over a heterogeneous family. Raw channel maps should still be returned for interpretation.

Implementation detail: the first version can store null channel arrays for moderate tests. Add streaming online mean/variance plus max-null accumulation before scaling to large production runs.

Production implementation should be explicitly two-pass:

1. Generate deterministic null channel maps to estimate `mu_null(rho, v, c)` and `sd_null(rho, v, c)`.
2. Regenerate the same deterministic null maps and compute the standardized max-null distribution.

This avoids storing huge null arrays while preserving studentization.

## FWER Correction

For every permutation, compute the maximum over all requested scales, locations, and channels:

```text
M_b = max_{rho, v, c} T_b(rho, v, c)
```

Corrected p-values use the package's plus-one convention:

```text
p_FWER = (1 + sum(M_b >= T_obs)) / (B + 1)
```

Report:

- `omnibus_p`: minimum corrected p-value over the whole family.
- `best_scale`: scale with minimum corrected p-value.
- `finest_detected_scale`: smallest FWHM with any corrected point at `alpha`.
- `stable_supported_range`: neighboring supported scales or a connected scale-space ridge.
- `supported`: table by scale and channel.

Be explicit that a corrected point at `(rho, v)` means significant evidence for a localized test centered at `v` with bandwidth `rho`, not voxelwise significance.

Do not assume scale support is monotone. Report an isolated finest detected scale, but prefer stable supported ranges for interpretation when neighboring scales agree.

## Result Object

Return class `scale_space_energy_result`:

```r
list(
  method = "scale_space_energy",
  observed = list(raw = ..., statistic = ...),
  p_corr = list(...),
  max_null = numeric(n_perm),
  scales = scale_space_spec,
  channels = character(),
  summary = data.frame(...),
  selected = list(best = ..., finest_detected = ..., stable_range = ...),
  mask = mask,
  params = list(...)
)
```

Add:

- `print.scale_space_energy_result()`
- `summary.scale_space_energy_result()`
- `plot_resolution_curve()`
- `plot_scale_space_map()`
- `energy_mass_map()`

Keep energy-mass maps descriptive unless they are derived from corrected scale-space hits.

For local energy, return three map types:

- `raw_energy = S_rho z^2`
- `excess_energy = S_rho z^2 - mu_null(rho, v)`
- `studentized = [S_rho z^2 - mu_null(rho, v)] / sd_null(rho, v)`

Energy-mass summaries should use above-null maps, such as `excess_energy` or the positive part of the studentized field, not raw `S_rho z^2`.

## Implementation Phases

### Phase 1: Internal Contracts

- Add image matrix coercion and mask/coordinate extraction.
- Add deterministic one-sample and two-sample statistic-map functions.
- Add within-subject contrast-image construction for repeated-measures designs.
- Add deterministic resampling index/sign generators.
- Add tests for shapes, finite output, mask consistency, and reproducibility.

### Phase 2: Scale-Space Core

- Add `scale_space()` and smoother construction.
- Implement `grid_gaussian` first.
- Add identity-scale support for voxelwise endpoint.
- Add effective-basis estimates, initially exact for small cases and approximate for larger cases.

### Phase 3: Low-Level Inference

- Implement `scale_space_energy(stat_vol, null_fun, ...)`.
- Compute raw and studentized observed channels.
- Implement two-pass studentization and max-null accumulation.
- Add corrected p-value maps.
- Add summary table with best scale, finest detected scale, and stable supported range.

### Phase 4: Raw-Image Wrappers

- Add `scale_space_onesample()`.
- Add `scale_space_twosample()`.
- Add `make_within_subject_contrast_images()`.
- Keep both wrappers thin: they should build `stat_vol`, `mask`, and `null_fun`, then call the low-level scanner.

### Phase 5: GLM Resampling

- Add `make_null_fun_glm()` with Freedman-Lane.
- Support formulas, numeric design matrices, named contrasts, and contrast vectors.
- Add two-group-with-covariates through this path rather than duplicating logic.
- Document this as the preferred route for unbalanced clinical/control comparisons.
- Add roadmap notes for heteroskedasticity-robust residual resampling.

### Phase 6: Graph Backend

- Add `graph_heat` using sparse mask graph construction.
- Borrow/adapt `adjoin` adjacency and Laplacian code with minimal dependencies.
- Add tests against dense `expm` on tiny graphs.
- Add optional band-pass/wavelet channel only after the low-pass scale-space path is stable.

### Phase 7: Reporting and Vignettes

- Add resolution-evidence curve.
- Add selected-scale maps and energy-mass summaries.
- Add vignette comparing global energy, scale-space localized energy, TFCE, LISA-style post-statistic filtering precedent, and voxelwise max-T on simulated same-sign and mixed-sign patterns.

## Test Plan

Follow contract-first computational testing.

### Contract Tests

- `as_image_matrix()` preserves subject count, mask length, dimensions, and template.
- `null_fun(b)` always returns length `n_mask`.
- Repeated calls to the same permutation index return identical maps.
- Invalid group labels, mismatched image dimensions, empty masks, and nonfinite covariates fail early.

### Differential Tests

- One-sample statistic maps match `t.test()` voxel by voxel on tiny fixtures.
- Two-sample statistic maps match base R calculations for Welch and pooled t.
- Tiny graph heat smoothing matches dense matrix exponential.
- Identity scale matches unsmoothed channels exactly.
- Global endpoint matches the chosen global statistic, initially `sum(z^2)` for local energy.

### Metamorphic Tests

- Reordering subjects leaves observed maps unchanged and permuted maps equivalent after reindexing.
- Multiplying all one-sample images by `-1` flips signed coherence but leaves local energy unchanged.
- Increasing FWHM should not increase effective basis count.
- Adding an all-zero subject effect leaves finite maps and stable mask handling.
- Repeated-measures contrast images match hand-computed cell contrasts on tiny fixtures.

### Null Calibration Tests

- One-sample sign-flip null simulations are not anti-conservative under symmetric nulls.
- Two-group label permutation is not anti-conservative under exchangeable equal-null groups.
- Plus-one corrected p-values and thresholds match `fwer_p_from_maxnull()` behavior.

### Power/Behavior Tests

- Same-sign diffuse blobs are detected by signed coherence at coarse scales.
- Mixed-sign distributed blobs are detected by local energy when signed coherence weakens.
- Voxel-focal blobs prefer small scales or identity endpoint.

### Performance Guardrails

- Record baseline runtime for small, medium, and representative masks.
- Track memory for null storage versus streaming mode.
- Keep dense exact backends test-only or tiny-mask-only.

## Open Design Decisions

- Whether `grid_gaussian` should be the default even for irregular masks, with graph heat as an opt-in.
- Whether studentization should be required for multi-channel scans or only the default. The working assumption is that it is required for inference and optional only for descriptive maps.
- How to estimate effective basis count cheaply for large masks and non-symmetric smoothers.
- How to define stable supported ranges: adjacent supported scales, connected scale-space ridges, or both.
- Whether `positive_energy` belongs in the default family or should require explicit prespecification because it adds multiplicity. The working assumption is explicit prespecification.

## First Milestone Definition

The first reviewable milestone is complete when:

- `scale_space_energy()` runs from `stat_vol + null_fun`.
- identity endpoint, global endpoint, and Gaussian-grid scales work.
- signed coherence and local energy channels work.
- two-pass studentization is implemented.
- corrected p-values control over scale, location, and channel.
- `summary()` reports best scale, finest detected scale, stable supported range, global endpoint, and voxelwise endpoint.
- `plot_resolution_curve()` works.
- tests cover contracts, exact small-case statistics, identity/global endpoints, plus-one FWER, and basic null calibration.

Raw-image wrappers, repeated-measures contrast helpers, GLM/Freedman-Lane support, graph heat, positive energy, and band-pass channels follow in later milestones.

## Suggested Reporting Language

We tested planned contrast images using scale-space localized energy inference. For each observed and permuted statistic map, we computed signed-coherence and local-energy fields over a prespecified Gaussian scale-space, including a global-energy endpoint and an identity voxelwise endpoint. Familywise error was controlled over all spatial scales, locations, and statistic channels using the permutation maximum. We report the scale with strongest corrected evidence, the finest spatial resolution at which the effect was detected, and the stable supported scale range. A significant result at scale `rho` indicates evidence for interaction energy localized to a neighborhood of bandwidth `rho`; it is not a voxelwise claim for every voxel in that neighborhood.

## LISA Precedent

LISA is useful precedent because it performs post-statistic spatial filtering of z-maps and applies the same transformation to permutation maps before inference. It also includes a generic mode where users supply an observed map and permuted maps. The distinction is important: LISA uses an adaptive local filter and FDR scores in filtered maps; scale-space localized energy inference searches a prespecified resolution family and uses max-null FWER across scale, location, and channel.

Reference: Lohmann et al. (2018), "LISA improves statistical analysis for fMRI", Nature Communications, https://www.nature.com/articles/s41467-018-06304-z.

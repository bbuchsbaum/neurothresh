# API Consistency Plan

Status: proposal — no backward compatibility required.
Goal: every public inference function in `neurothresh` follows one regular shape, with one null contract, one tail vocabulary, one stat-map name, and one return type.

## Principles

1. **One name per concept.** No `z_vol` vs `stat_vol`, no `two_sided=TRUE` vs `tail="two"`, no `fwhm_vox` vs `fwhm_mm` competing in the same signature.
2. **Pluggable null is the universal seam.** Every permutation-based method takes the same `null_spec` argument. Parametric methods (RFT) take `smoothness` instead.
3. **Caller owns canonicalization.** Inference functions assume the input is already a Z-equivalent map. `canonicalize()` is a separate exported step.
4. **One result class.** Every public method returns an S3 `neurothresh_result` with predictable fields.

## 1. Canonical signature shape

```r
method_name(
  stat_vol,                       # 3D Z-equivalent map (always this name)
  mask        = NULL,             # logical 3D, defaults to is.finite(stat_vol)
  alpha       = 0.05,             # (or `q` for FDR methods only)
  tail        = c("pos", "neg", "two"),
  null_spec   = NULL,             # permutation methods only
  n_perm      = 5000,             # permutation methods only
  smoothness  = NULL,             # RFT / mc-based methods only
  seed        = NULL,
  ...                             # method-specific knobs (kappa, H/E, gamma, ...)
)
```

Drop forever: `z_vol`, `two_sided`, `two_sided_policy`, `fwhm_vox`, `fwhm_mm` (as a top-level arg), `acf_params`, `stat_type`, `df`, `allow_voxel_signflip`, `null = c("signflip_voxel", ...)` string sugar.

## 2. The `null_spec` contract

A `null_spec` is an S3 object (`class = "neurothresh_null_spec"`) — a list with:

| Field      | Type                          | Meaning                                                 |
|------------|-------------------------------|---------------------------------------------------------|
| `stat_vol` | 3D array / NeuroVol           | Observed Z-equivalent map.                              |
| `mask`    | 3D logical                    | Mask used to build the null.                            |
| `draw`    | `function(b) -> numeric(n_mask)` | Returns one null draw, mask-ordered, 1 ≤ b ≤ n_perm. |
| `n_perm`  | integer                       | Maximum supported draws.                                |
| `df`      | numeric (optional)            | Degrees of freedom (subject-level nulls).               |
| `kind`    | character                     | `"subject_signflip"`, `"two_sample_perm"`, `"voxel_signflip"`, `"mc_fwhm"`, `"mc_acf"`, `"parcel_perm"`. |
| `meta`    | list                          | Free-form (smoothness used, etc.).                      |

All factories return one of these. Today's divergence — `make_null_fun_voxel_signflip` / `make_null_fun_mc_fwhm` / `make_null_fun_mc_acf` return bare functions, while `make_null_fun_subject_signflip` / `make_null_fun_onesample` / `make_null_fun_twosample` return rich lists — collapses to the single shape above. The `$draw` slot replaces the current mix of `$null_fun` / bare function returns.

### Rename factory family

```
null_subject_signflip(subject_vols, mask = NULL, n_perm = 5000, seed = NULL, ...)
null_two_sample(subject_vols, group, mask = NULL, n_perm = 5000, seed = NULL, ...)
null_voxel_signflip(stat_vol, mask = NULL, n_perm = 5000, seed = NULL)
null_mc_fwhm(template, mask, fwhm_mm, n_perm = 5000, seed = NULL)
null_mc_acf(template, mask, acf, n_perm = 5000, seed = NULL)
null_parcel_perm(stat_vol, parcels, mask = NULL, n_perm = 5000, seed = NULL)   # new
```

Drop the `make_null_fun_*` prefix; verb-noun is shorter and `null_*` reads well at the call site:

```r
ns  <- null_subject_signflip(subject_vols, n_perm = 5000)
res <- octree_scan_stepdown(ns$stat_vol, mask = ns$mask, null_spec = ns, alpha = 0.05)
```

Inference functions accept `null_spec` and pull `stat_vol`/`mask` from it when not separately supplied. If the caller passes `stat_vol` directly without a `null_spec`, the function errors with a clear message pointing at the factory list (no silent fallback to voxel sign-flip).

## 3. The `smoothness` contract (parametric / MC paths)

Single tagged-list constructor family, one canonical unit (mm):

```r
smoothness_fwhm(fwhm_mm)                # length-1 or length-3 numeric, mm
smoothness_acf(a, b, c)                 # AFNI-style ACF parameters
```

Both return an object with `class = "neurothresh_smoothness"` and known slots. RFT methods accept only `smoothness_fwhm()`. MC null factories accept either. The package exposes no other smoothness arguments.

## 4. The `tail` contract

`tail = c("pos", "neg", "two")` everywhere. `"two"` always means: combine via per-voxel max over positive and negative tails (TFCE) or symmetric two-sided p-values (RFT, cluster-FDR, scans). Drop `two_sided_policy` from `cluster_fdr`; the `"two"` semantics for FDR collapses to the previously-named `BH_all` policy. (The `split_q` variant, if anyone wants it, becomes a separate function `cluster_fdr_split` rather than an argument.)

## 5. Canonicalization is a separate step

```r
canonicalize(stat_vol, stat_type = c("Z", "t", "log10p"), df = NULL, ...) -> stat_vol
```

`hier_scan` loses `stat_type` and `df`. Users who supply a t-map call `canonicalize()` first. Documentation and vignettes update to lead with this two-line idiom:

```r
zvol <- canonicalize(t_map, stat_type = "t", df = 24)
res  <- hier_scan(zvol, parcels = atlas, alpha = 0.05)
```

## 6. The result class

Every public method returns `class = "neurothresh_result"` with consistent slots:

| Field        | Always present? | Meaning                                                        |
|--------------|-----------------|----------------------------------------------------------------|
| `stat_vol`   | yes             | Observed map (canonicalized).                                  |
| `mask`       | yes             | Logical mask used.                                             |
| `method`     | yes             | One of `"hier_scan"`, `"octree_scan_fwer"`, `"octree_scan_stepdown"`, `"tfce_fwer"`, `"cluster_fdr"`, `"rft_peak_fwer"`, `"rft_cluster_fwer"`. |
| `threshold`  | yes             | Critical value on the relevant statistic.                      |
| `reject`     | yes             | 3D logical volume of rejected voxels.                          |
| `pvals`      | yes             | 3D corrected p-value volume.                                   |
| `clusters`   | when applicable | Data frame of clusters (size, peak, p, label).                 |
| `null_max`   | permutation only| Vector of max-null statistics.                                 |
| `params`     | yes             | List capturing alpha, tail, n_perm, kappa, H/E, smoothness, etc. |
| `call`       | yes             | `match.call()`.                                                |

Print/summary/plot methods key off `method` and these uniform slots.

## 7. Per-function signature targets

```r
hier_scan(stat_vol, parcels, prior_vol = NULL, mask = NULL,
          alpha = 0.05, tail = c("pos","neg","two"),
          kappa = c(0.5, 1, 2), prior_eta = 0.9,
          null_spec = NULL, n_perm = 5000,
          method = c("stepdown","hierarchical"),
          gamma = 0.5, gamma_root = 0.5, min_voxels = 8, min_alpha = 1e-6,
          parallel = FALSE, n_workers = NULL, seed = NULL)

octree_scan_fwer(stat_vol, prior_vol = NULL, mask = NULL,
                 alpha = 0.05, tail = c("pos","neg","two"),
                 prior_eta = 0.9,
                 null_spec = NULL, n_perm = 5000, seed = NULL)

octree_scan_stepdown(stat_vol, prior_vol = NULL, mask = NULL,
                     alpha = 0.05, tail = c("pos","neg","two"),
                     prior_eta = 0.9, gamma = 0.5,
                     min_voxels = 8, min_alpha = 1e-6,
                     report = c("coarsest","finest","greedy"),
                     null_spec = NULL, n_perm = 5000, seed = NULL)

tfce_fwer(stat_vol, mask = NULL,
          alpha = 0.05, tail = c("pos","neg","two"),
          H = 2.0, E = 0.5, dh = 0.1,
          null_spec = NULL, n_perm = 5000, seed = NULL)

cluster_fdr(stat_vol, mask = NULL,
            q = 0.05, tail = c("pos","neg","two"),
            cluster_thresh = 3.0,
            p_method = c("rft","perm"),
            smoothness = NULL,                  # required if p_method = "rft"
            null_spec = NULL,                   # required if p_method = "perm"
            n_perm = 5000, seed = NULL)

rft_peak_fwer(stat_vol, mask = NULL,
              alpha = 0.05, tail = c("pos","neg","two"),
              smoothness, df = Inf)             # smoothness via smoothness_fwhm()

rft_cluster_fwer(stat_vol, mask = NULL,
                 alpha = 0.05, tail = c("pos","neg","two"),
                 cluster_thresh = 3.0,
                 smoothness, df = Inf)

tfce_transform(stat_vol, mask = NULL,
               H = 2.0, E = 0.5, dh = 0.1,
               tail = c("pos","neg","two"))    # gain "two"
```

`rft_fwer_twosided` disappears; `tail = "two"` on `rft_peak_fwer` / `rft_cluster_fwer` covers it.

## 8. Implementation phases

1. **Spec + null_spec class.** Create `R/null_spec.R`: define class, validators, conversion helpers. Update `make_null_fun_*` factories to return `null_spec` objects (rename to `null_*`).
2. **Smoothness class.** Create `R/smoothness.R`: `smoothness_fwhm()`, `smoothness_acf()`, validators.
3. **Canonicalize export.** Promote `canonicalize()` (internal in hier_scan today) to an exported function; move `stat_type`/`df` logic out of `hier_scan`.
4. **Migrate inference functions one at a time** — `octree_scan_fwer` → `octree_scan_stepdown` → `tfce_fwer` → `cluster_fdr` → `hier_scan` → `rft_*`. Each migration: update signature, drop renamed args, route through `null_spec` / `smoothness`, return `neurothresh_result`.
5. **Result class.** Add `R/result.R` with constructor, `print.neurothresh_result`, `summary.neurothresh_result`, `as.data.frame.neurothresh_result` for clusters. Replace ad-hoc returns.
6. **Tests.** Rewrite tests against new signatures. Keep numerical fixtures (seeds, expected p-values) — only the calling convention changes.
7. **Docs + vignettes.** Update every `@examples` block, `vignettes/neurothresh-intro.Rmd`, and `vignettes/simulation-validation.Rmd` to the new shape. Update `AGENTS.md` / `CLAUDE.md` algorithm-flow snippets if they reference old arg names.
8. **NEWS.** Single big "0.x.0 — API regularization" entry listing every renamed function and argument with a one-line porting recipe.

## 9. Polymorphic input layer (Phase 1b)

The package already has the seed of this layer as `as_image_matrix()` in `R/scale_space_utils.R:15`, used by `null_one_sample` / `null_two_sample` / `null_glm` factories. Rather than introduce a parallel abstraction, **generalize and rename it** to `as_subject_matrix()`, broaden the input contract to cover neuroim2 4D types, and route every subject-level null factory through it.

### 9.1 Function

```r
as_subject_matrix(x, mask = NULL, template = NULL, ...) -> structure(list(
  X        = matrix,         # n_subjects × n_mask, finite, double
  mask     = logical array,  # 3D logical, possibly LogicalNeuroVol if input carried space
  mask_idx = integer,        # which(mask), 1-based
  template = array/NeuroVol, # 3D template (NeuroVol when input had spatial metadata)
  dims     = integer(3),
  coords   = matrix,         # n_mask × 3 array indices
  space    = NULL | NeuroSpace,
  source   = character       # "list", "neurovec", "sparseneurovec", "matrix"
), class = "neurothresh_image_matrix")
```

Promises are typed only as strongly as the input supports: `mask`/`template` upgrade to `LogicalNeuroVol`/`NeuroVol` when the input carries `NeuroSpace`, otherwise they stay plain arrays. `space` is `NULL` for plain-array inputs.

### 9.2 Input branches

Implement as a plain exported function with type branches (no S4 in the first slice; if internal complexity grows, split dispatch behind methods later):

| Branch                         | Behavior                                                                                         |
|--------------------------------|--------------------------------------------------------------------------------------------------|
| `list<array>` / `list<NeuroVol>` | Stack as rows; require identical dims; preserve `NeuroSpace` from `images[[1]]` if present.    |
| `NeuroVec`                     | Subjects on 4th dim. **Note**: `as.matrix(NeuroVec)` is `n_voxels × n_volumes` — must transpose, then mask-subset, to get `n_subjects × n_mask`. |
| `SparseNeuroVec`               | Use the intrinsic mask; supplied `mask` may further restrict (see precedence below).             |
| `matrix`                       | Already `n_subjects × n_voxels`; `mask` and/or `template` must be supplied.                      |

### 9.3 Mask precedence (strict)

- `mask = NULL`: infer the finite intersection across subjects. For `SparseNeuroVec`, use its intrinsic mask.
- `mask` supplied: **use exactly that domain**. Error if it includes non-finite values or voxels outside a `SparseNeuroVec`'s intrinsic mask. No silent intersection. (A future opt-in `strict_mask = FALSE` may relax this — out of scope here.)

### 9.4 Null factory naming (design-facing)

Statistical design names the factory; mechanism is an implementation detail:

```r
null_one_sample(images, mask = NULL, n_perm = 5000, seed = NULL, z_equiv = TRUE)
null_two_sample(images, group, mask = NULL, levels = NULL, contrast = NULL,
                variance = c("welch", "pooled"), n_perm = 5000, seed = NULL, z_equiv = TRUE)
null_glm(images, design, contrast, mask = NULL, n_perm = 5000, seed = NULL, z_equiv = TRUE)
null_subject_signflip(...)   # alias of null_one_sample for users who want the mechanism name
```

Two-sample contrast direction: default is `level 2 minus level 1`. Documented explicitly. Users can pin it with `levels = c("control", "case")` or `contrast = c("case", "control")`.

### 9.5 Consolidation of duplicate factories

`make_null_fun_subject_signflip()` (samples signs lazily per call, no `n_perm`) and `make_null_fun_onesample()` (precomputes a deterministic `n_perm × n_subj` sign matrix indexed by `b`) are near-duplicates. The deterministic-by-index contract is correct for `null_spec$draw(b)` (reproducible permutation indexing, supports caching, supports parallel resampling). **Keep that one.** The lazy variant is retired; `null_subject_signflip` becomes a thin alias of `null_one_sample`.

### 9.6 Acceptance test

Single canonical test, run for each input type:

```r
toy_array_list <- ...                          # list of 3D arrays
toy_neurovols  <- lapply(toy_array_list, as_neurovol)
toy_neurovec   <- as_neurovec(toy_array_list)
toy_sparse     <- as_sparseneurovec(toy_array_list, mask = toy_mask)
toy_matrix     <- subject_voxel_matrix(...)

a <- as_subject_matrix(toy_array_list, mask = toy_mask)
b <- as_subject_matrix(toy_neurovols,  mask = toy_mask)
c <- as_subject_matrix(toy_neurovec,   mask = toy_mask)
d <- as_subject_matrix(toy_sparse)             # mask intrinsic
e <- as_subject_matrix(toy_matrix, mask = toy_mask, template = toy_template)

stopifnot(identical(a$X, b$X), identical(a$X, c$X),
          identical(a$X, d$X), identical(a$X, e$X))
stopifnot(identical(a$mask_idx, b$mask_idx), ..., identical(a$mask_idx, e$mask_idx))
# Reconstructed stat_vol from null_one_sample(.)$stat_vol must match across inputs.
```

### 9.7 Dependency relationship

Phase 1b is a **prerequisite** for Phase 1 factory rename and for every Phase 4 migration that builds nulls from subject data. It is independent of Phase 2 (smoothness) and Phase 5 (result class) and can be done in parallel with them.

## 10. Out of scope (intentionally)

- Renaming the methods themselves (`hier_scan`, `octree_scan_*`, `tfce_fwer`, `cluster_fdr`, `rft_*`) — they're already verb-y and method-specific.
- Changing the underlying statistics or any numerical behavior.
- Touching Rcpp signatures — only the R-facing layer changes.
- Adding new methods. This is purely shape work.

## 11. Acceptance criteria

- `grep -nE "z_vol|two_sided|two_sided_policy|fwhm_vox|allow_voxel_signflip|stat_type" R/` returns nothing in the public API surface.
- Every exported inference function matches the canonical signature shape from §1.
- Every permutation-based method accepts (and only accepts) `null_spec`.
- `R CMD check --as-cran .` is clean.
- A single short vignette demonstrates: build a `null_spec`, build a `smoothness` (where needed), run each method, print/summary the result — same idiom for all of them.

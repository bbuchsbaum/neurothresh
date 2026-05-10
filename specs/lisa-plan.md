# LISA Local Adaptive Inference Plan

## Positioning

Implement LISA as a separate local-adaptive inference engine, not as a
scale-space energy channel.

- `scale_space_energy()` asks: what spatial resolution supports a distributed
  effect?
- `lisa()` asks: which voxels have locally supported evidence after adaptive
  post-statistic filtering?

Both methods should share the same upstream map/resampling infrastructure, but
they make different inferential claims. LISA controls FDR in the filtered
statistic space; it does not provide FWER, voxelwise raw-statistic significance,
or a finest-supported spatial scale.

## Refactored Package Contract

The current implementation direction makes `neurothresh_null_spec` the universal
permutation seam. LISA should therefore use:

```r
fit <- lisa(
  null_spec = ns,
  filter = lisa_bilateral(),
  tail = "two",
  q = 0.05
)
```

The default path pulls `stat_vol`, `mask`, `draw`, and `n_perm` from
`null_spec`. Direct `stat_vol`/`mask` arguments may be allowed as explicit
overrides, but a calibrated null source is required.

Inputs are assumed to be canonicalized Z-equivalent maps. If callers have
t-statistics or `-log10(p)` maps, they should call `canonicalize()` first or use
a null factory that already produces Z-equivalent maps.

Use package-wide conventions:

- `tail = c("pos", "neg", "two")`
- raw-image wrappers call `null_one_sample()`, `null_two_sample()`, or
  `null_glm()`
- external observed/null maps enter through a `null_spec` factory
- public inference returns `neurothresh_result`

## Public API

Core engine:

```r
lisa <- function(null_spec,
                 stat_vol = NULL,
                 mask = NULL,
                 filter = lisa_bilateral(),
                 tail = c("pos", "neg", "two"),
                 q = 0.05,
                 p0 = 1,
                 monotone = TRUE,
                 scale_n_null = 30,
                 backend = c("neuroim2", "reference", "adjoin_sparse"))
```

Filter specification:

```r
lisa_bilateral(
  radius = 2,
  sigma_space = 2,
  sigma_range = 2,
  n_iter = 2,
  scale = "sd_first_30_nulls",
  boundary = c("neuroim2", "mask_normalized", "lisa2018")
)
```

`neurothresh` should not own a new general bilateral smoother. The first
production backend should adapt `neuroim2::bilateral_filter()`, because
`neuroim2` is already an imported package and provides compiled 3D bilateral
filtering. The `lisa_bilateral()` object owns LISA-specific choices: filter
defaults, range scaling from observed/null maps, repeated application over null
draws, boundary policy labeling, and diagnostics.

The `adjoin_sparse` backend is an optional graph formulation. It can use
`adjoin::bilateral_smoother()`-style sparse weights over masked voxel
coordinates and statistic features. This is elegant and mask-native, but the
weights depend on the statistic map, so they must be rebuilt for each observed
and null map. Treat it as a future/backend comparison path, not the first
required implementation.

External null-map bridge:

```r
null_from_maps(stat_vol, null_maps, mask = NULL,
               kind = "external_maps", meta = list())
```

Wrappers:

```r
lisa_onesample(images, mask = NULL, n_perm = 5000, ...)
lisa_twosample(images, group, mask = NULL, n_perm = 5000, ...)
lisa_glm(images, design, contrast, mask = NULL, n_perm = 5000, ...)
```

Wrappers should construct a `null_spec` and then call `lisa()`. They should not
own filtering, FDR estimation, or result construction.

## Result Object

LISA should return `neurothresh_result`:

- `method = "lisa"`
- `stat_vol`: observed canonicalized map
- `mask`: analysis mask
- `threshold`: filtered-statistic cutoff at requested `q`
- `reject`: logical discovery map
- `pvals`: FDR/q-value map in filtered-statistic space
- `params`: `q`, `tail`, `n_perm`, `p0`, filter settings, backend
- `extras`: filtered statistic map, FDR curve, histograms, scaling diagnostics

The filtered statistic should be stored in `extras$filtered_stat`; it is not a
standardized Z-map after bilateral filtering.

## Algorithm

1. Validate `null_spec` and resolve `stat_vol`, `mask`, `draw`, and `n_perm`.
2. Estimate the bilateral range scaling from the observed map and/or first
   `scale_n_null` null draws, matching the LISA default of using early null-map
   variability for range calibration.
3. Apply the selected bilateral backend to the observed statistic map.
4. Apply the same backend and scaling policy to each null draw, streaming where
   possible.
5. Estimate empirical filtered-null tail probabilities from all filtered null
   voxel values.
6. Estimate empirical mixture tail probabilities from the filtered observed map.
7. Compute conservative local FDR/q maps with `p0 = 1` by default.
8. Apply monotonicity correction so stronger filtered evidence cannot receive
   worse FDR values than weaker evidence in the same tail.
9. Build `neurothresh_result`.

Tail handling is channel-specific:

- `pos`: use `Lambda >= lambda`
- `neg`: use `Lambda <= lambda`
- `two`: use `abs(Lambda) >= abs(lambda)`

## Implementation Layers

### Layer 1: Filter Backend Adapter

Add `R/lisa_filter.R`:

- `lisa_bilateral()` filter spec constructor
- `is_lisa_filter()`
- `lisa_apply_filter()` internal dispatcher
- `.lisa_filter_neuroim2()` adapter around `neuroim2::bilateral_filter()`
- `.lisa_filter_reference()` small audited fallback for fixtures and boundary
  tests
- `.lisa_filter_adjoin_sparse()` optional graph backend if `adjoin` is
  available

The first version should prioritize an adapter around the existing compiled
`neuroim2` implementation, not a new C++ smoother in `neurothresh`. The
reference path should exist only to define expected behavior on tiny arrays and
to test mask/boundary semantics.

Boundary behavior must be explicit. `neuroim2::bilateral_filter()` filters
masked output voxels but reads neighborhood values from the full array; outside
mask voxels can therefore affect boundary voxels depending on how the input
array is filled. Start with `boundary = "neuroim2"` as the faithful adapter
mode, add `boundary = "mask_normalized"` to exclude out-of-mask neighbors, and
reserve `boundary = "lisa2018"` for the paper-style median fallback/discard
rules.

### Layer 2: FDR Backend

Add `R/lisa_fdr.R`:

- `lisa_fdr_map(filtered_obs, filtered_nulls, mask, tail, q, p0, monotone)`
- internal monotone empirical FDR helpers
- histogram/streaming accumulator once correctness is locked

The first implementation may store filtered null vectors for small tests, but
the production target is streaming null histograms to avoid `B x V` memory use.

### Layer 3: Inference Engine

Add `R/lisa.R`:

- `lisa()`
- result construction via `new_neurothresh_result()`
- `filtered_stat()` or `lisa_filtered_stat()` extractor if needed
- `plot_lisa_map()` only if it follows existing map-extraction conventions

### Layer 4: Wrappers

Add `R/lisa_wrappers.R`:

- `lisa_onesample()` -> `null_one_sample()` -> `lisa()`
- `lisa_twosample()` -> `null_two_sample()` -> `lisa()`
- `lisa_glm()` -> `null_glm()` -> `lisa()`

For unbalanced clinical/control designs, documentation should recommend
`lisa_glm()` over simple two-sample label permutation.

## Validation Plan

Mechanical tests:

- constant maps remain constant after filtering
- impulse maps are locally smoothed
- mask boundaries do not leak outside the mask
- `backend = "neuroim2"` matches direct `neuroim2::bilateral_filter()` under
  `boundary = "neuroim2"`
- reference filter output is deterministic on tiny arrays
- `backend = "adjoin_sparse"` is optional and skips cleanly when `adjoin` is not
  installed
- `null_from_maps()` validates null-map dimensions and draw lengths
- `lisa()` requires `null_spec`
- `tail = pos/neg/two` produce expected monotone ordering
- FDR maps are bounded in `[0, 1]`
- `threshold` and `reject` agree with `q`
- result has class `neurothresh_result`

Statistical tests:

- all-null simulations have conservative/acceptable empirical FDR
- local coherent signal favors LISA over max-T in sensitivity examples
- mixed-sign distributed signal is not overclaimed by LISA
- same `null_spec` can feed both `lisa()` and `scale_space_energy()`

## Documentation

Add a vignette or section named "LISA local adaptive inference":

- LISA is post-statistic adaptive filtering plus FDR in filtered space.
- Filtered values are not Z-scores.
- LISA discoveries are locally supported filtered-map discoveries, not raw
  voxelwise discoveries.
- `scale_space_energy()` and `lisa()` are complementary and should often be run
  from the same `null_spec`.

Example:

```r
ns <- null_glm(
  images = D_task_model,
  design = ~ group + age + sex + motion,
  contrast = "group",
  mask = gm_mask,
  n_perm = 10000,
  statistic = "z"
)

fit_sse <- scale_space_energy(null_spec = ns)
fit_lisa <- lisa(null_spec = ns, tail = "two", q = 0.05)
```

The `scale_space_energy(null_spec = ns)` call requires migrating the scale-space
engine to the same `null_spec` interface or adding a compatibility adapter.

## Backend Decision

LISA needs bilateral filtering to be faithful to the method, but `neurothresh`
does not need to implement a full bilateral image smoother. The package should
own the LISA inference contract and delegate the filtering primitive:

1. First backend: `neuroim2::bilateral_filter()`.
2. Test/reference backend: small R implementation for audited fixtures.
3. Optional backend: sparse graph smoother borrowed from `adjoin`.

This keeps the package focused on inference while reusing existing image and
graph infrastructure.

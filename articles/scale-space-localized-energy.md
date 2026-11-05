# Scale-Space Localized Energy Inference

[`scale_space_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md)
tests a planned statistic map across a prespecified resolution
continuum. The low-level API deliberately follows the rest of
`neurothresh`: the inferential engine consumes an observed statistic map
plus a calibrated null generator.

There are two related but distinct questions:

- **Localization:** at what resolution can an effect be detected?
- **Frequency readout:** at what spatial band is the effect strongest?

[`local_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/local_energy.md)
answers the first question. For the second, use
[`global_lowpass_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/global_lowpass_energy.md)
or
[`global_band_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/global_band_energy.md).

``` r

library(neurothresh)

set.seed(1)
dims <- c(12, 12, 12)
signal_mask <- array(FALSE, dims)
signal_mask[4:9, 4:9, 4:9] <- TRUE
subject_images <- replicate(
  8,
  {
    img <- array(rnorm(prod(dims), sd = 0.5), dim = dims)
    img[signal_mask] <- img[signal_mask] + 3
    img
  },
  simplify = FALSE
)
```

For a one-sample planned contrast, the raw-image wrapper builds the
group map and deterministic subject-level sign-flip null maps, then
calls the same low-level scanner. The `n_perm` values used below are
kept tiny so the vignette builds quickly: with `B = 19` permutations the
smallest possible corrected p-value is `1/(B + 1) = 0.05`, which is far
too coarse for real analyses. Production runs should use at least
`n_perm = 1000`, and 5,000 to 10,000 if power is a concern.

``` r

fit <- scale_space_onesample(
  subject_images,
  scales = scale_space(
    fwhm_mm = c(4, 8, 12),
    include_global = TRUE,
    include_identity = TRUE
  ),
  statistics = list(
    global_lowpass_energy(),
    global_band_energy(),
    local_energy(),
    signed_coherence(tail = "two")
  ),
  n_perm = 19,
  seed = 1,
  alpha = 0.1
)

fit
#> Scale-space localized energy inference
#>   permutations: 19 
#>   channels: global_lowpass_energy, global_band_energy, local_energy, signed_coherence_two 
#>   omnibus p_FWER: 0.05 
#>   frequency-readout channels: global_lowpass_energy, global_band_energy
head(summary(fit))
#>   scale_index channel_index scale_type fwhm_mm               channel
#> 1           1             1     global     Inf global_lowpass_energy
#> 2           1             2     global     Inf    global_band_energy
#> 3           1             3     global     Inf          local_energy
#> 4           1             4     global     Inf  signed_coherence_two
#> 5           2             1   gaussian       4 global_lowpass_energy
#> 6           2             2   gaussian       4    global_band_energy
#>                           channel_role channel_scope frequency_readout min_p
#> 1  cumulative low-pass spectral energy        global              TRUE  0.05
#> 2     global band-pass spectral energy        global              TRUE  1.00
#> 3 localized squared-statistic evidence         local             FALSE  0.05
#> 4            signed low-pass coherence         local             FALSE  0.65
#> 5  cumulative low-pass spectral energy        global              TRUE  0.05
#> 6     global band-pass spectral energy        global              TRUE  0.05
#>   max_stat any_supported n_supported effective_basis p_fwer scale_label
#> 1 4.212548          TRUE           1               1   0.05      global
#> 2 0.000000         FALSE           0               1   1.00      global
#> 3 4.210267          TRUE           1               1   0.05      global
#> 4 3.868120         FALSE           0               1   0.65      global
#> 5 4.215862          TRUE           1              27   0.05        4 mm
#> 6 4.215932          TRUE           1              27   0.05        4 mm
#>                                                 curve_use
#> 1                               spatial-frequency readout
#> 2                               spatial-frequency readout
#> 3 localized squared evidence; not a spectral-energy curve
#> 4    signed low-pass coherence; opposite signs can cancel
#> 5                               spatial-frequency readout
#> 6                               spatial-frequency readout
```

The result controls familywise error over spatial scale, location, and
requested statistic channels. The global endpoint is a distributed
energy test. The identity endpoint is the voxelwise limit. Intermediate
FWHM values ask whether the effect is detectable at that localization
scale.

``` r

resolution <- plot_resolution_curve(fit)
resolution[, c("scale_label", "channel", "min_p", "any_supported", "curve_use")]
#>    scale_label               channel min_p any_supported
#> 1       global global_lowpass_energy  0.05          TRUE
#> 2       global    global_band_energy  1.00         FALSE
#> 3       global          local_energy  0.05          TRUE
#> 4       global  signed_coherence_two  0.65         FALSE
#> 5        12 mm global_lowpass_energy  0.05          TRUE
#> 6        12 mm    global_band_energy  0.05          TRUE
#> 7        12 mm          local_energy  0.05          TRUE
#> 8        12 mm  signed_coherence_two  0.35         FALSE
#> 9         8 mm global_lowpass_energy  0.05          TRUE
#> 10        8 mm    global_band_energy  0.05          TRUE
#> 11        8 mm          local_energy  0.05          TRUE
#> 12        8 mm  signed_coherence_two  0.20         FALSE
#> 13        4 mm global_lowpass_energy  0.05          TRUE
#> 14        4 mm    global_band_energy  0.05          TRUE
#> 15        4 mm          local_energy  0.05          TRUE
#> 16        4 mm  signed_coherence_two  0.25         FALSE
#> 17       voxel global_lowpass_energy  0.05          TRUE
#> 18       voxel    global_band_energy  1.00         FALSE
#> 19       voxel          local_energy  0.05          TRUE
#> 20       voxel  signed_coherence_two  0.45         FALSE
#>                                                  curve_use
#> 1                                spatial-frequency readout
#> 2                                spatial-frequency readout
#> 3  localized squared evidence; not a spectral-energy curve
#> 4     signed low-pass coherence; opposite signs can cancel
#> 5                                spatial-frequency readout
#> 6                                spatial-frequency readout
#> 7  localized squared evidence; not a spectral-energy curve
#> 8     signed low-pass coherence; opposite signs can cancel
#> 9                                spatial-frequency readout
#> 10                               spatial-frequency readout
#> 11 localized squared evidence; not a spectral-energy curve
#> 12    signed low-pass coherence; opposite signs can cancel
#> 13                               spatial-frequency readout
#> 14                               spatial-frequency readout
#> 15 localized squared evidence; not a spectral-energy curve
#> 16    signed low-pass coherence; opposite signs can cancel
#> 17                               spatial-frequency readout
#> 18                               spatial-frequency readout
#> 19 localized squared evidence; not a spectral-energy curve
#> 20    signed low-pass coherence; opposite signs can cancel
```

[`scale_space_diagnostics()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_diagnostics.md)
labels each curve by its intended interpretation. This is useful because
a local-energy curve should not be read as a pure spatial-frequency
spectrum.

``` r

scale_space_diagnostics(fit)[, c(
  "channel", "frequency_readout", "best_scale", "best_p_fwer",
  "supported_range", "curve_use"
)]
#>                 channel frequency_readout best_scale best_p_fwer
#> 1    global_band_energy              TRUE       4 mm        0.05
#> 2 global_lowpass_energy              TRUE     global        0.05
#> 3          local_energy             FALSE     global        0.05
#> 4  signed_coherence_two             FALSE       8 mm        0.20
#>   supported_range                                               curve_use
#> 1         4-12 mm                               spatial-frequency readout
#> 2         4-12 mm                               spatial-frequency readout
#> 3         4-12 mm localized squared evidence; not a spectral-energy curve
#> 4            <NA>    signed low-pass coherence; opposite signs can cancel
```

The selected scales are separated by target. `best_frequency` is the
strongest spatial-frequency readout, while `best_localization` is the
strongest location-wise scale-space evidence.

``` r

fit$selected$best_frequency[, c("channel", "scale_label", "p_fwer", "curve_use")]
#>                 channel scale_label p_fwer                 curve_use
#> 1 global_lowpass_energy      global   0.05 spatial-frequency readout
fit$selected$best_localization[, c("channel", "scale_label", "p_fwer", "curve_use")]
#>        channel scale_label p_fwer
#> 3 local_energy      global   0.05
#>                                                 curve_use
#> 3 localized squared evidence; not a spectral-energy curve
```

![Corrected evidence across the localization continuum. Higher values
indicate stronger familywise-corrected evidence; the dashed line marks
alpha =
0.10.](scale-space-localized-energy_files/figure-html/resolution-curve-plot-1.png)

Corrected evidence across the localization continuum. Higher values
indicate stronger familywise-corrected evidence; the dashed line marks
alpha = 0.10.

A corrected point at a scale does not mean each voxel is significant. It
means that a localized test centered at that location and bandwidth has
corrected evidence. For local-energy display, use above-null energy
summaries rather than raw `z^2` maps.

For spatial-frequency interpretation, prefer
[`global_band_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/global_band_energy.md):
its peak is the band with strongest corrected evidence.
[`global_lowpass_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/global_lowpass_energy.md)
is cumulative, so coarse curves can remain strong even when the effect
is not confined to one band.

``` r

mass <- energy_mass_map(fit)
range(as.numeric(mass), na.rm = TRUE)
#> [1] 0.0000000000 0.0006860179
```

![The statistic image contains a broad planted effect. The energy-mass
map shows where above-null scale-space energy is concentrated at the
selected
scale.](scale-space-localized-energy_files/figure-html/energy-mass-plot-1.png)

The statistic image contains a broad planted effect. The energy-mass map
shows where above-null scale-space energy is concentrated at the
selected scale.

For two-group analyses,
[`scale_space_twosample()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_twosample.md)
provides a simple label permutation wrapper. For unbalanced
clinical/control data or models with nuisance covariates, prefer a GLM
null generator.

``` r

group <- factor(rep(c("control", "case"), each = 4), levels = c("control", "case"))
age <- seq_along(group)

glm_null <- null_glm(
  subject_images,
  design = ~ group + age,
  data = data.frame(group = group, age = age),
  contrast = "groupcase",
  n_perm = 19,
  seed = 2
)

glm_fit <- scale_space_energy(
  stat_vol = glm_null$stat_vol,
  mask = glm_null$mask,
  null_spec = glm_null,
  scales = scale_space(c(4, 8), include_global = TRUE, include_identity = TRUE),
  statistics = list(global_band_energy(), local_energy()),
  n_perm = 19
)

summary(glm_fit)
#>   scale_index channel_index scale_type fwhm_mm            channel
#> 1           1             1     global     Inf global_band_energy
#> 2           1             2     global     Inf       local_energy
#> 3           2             1   gaussian       4 global_band_energy
#> 4           2             2   gaussian       4       local_energy
#> 5           3             1   gaussian       8 global_band_energy
#> 6           3             2   gaussian       8       local_energy
#> 7           4             1   identity       0 global_band_energy
#> 8           4             2   identity       0       local_energy
#>                           channel_role channel_scope frequency_readout min_p
#> 1     global band-pass spectral energy        global              TRUE   1.0
#> 2 localized squared-statistic evidence         local             FALSE   1.0
#> 3     global band-pass spectral energy        global              TRUE   1.0
#> 4 localized squared-statistic evidence         local             FALSE   1.0
#> 5     global band-pass spectral energy        global              TRUE   1.0
#> 6 localized squared-statistic evidence         local             FALSE   1.0
#> 7     global band-pass spectral energy        global              TRUE   1.0
#> 8 localized squared-statistic evidence         local             FALSE   0.8
#>      max_stat any_supported n_supported effective_basis p_fwer scale_label
#> 1  0.00000000         FALSE           0           1.000    1.0      global
#> 2 -0.30861189         FALSE           0           1.000    1.0      global
#> 3 -0.37095944         FALSE           0          27.000    1.0        4 mm
#> 4  0.94536592         FALSE           0          27.000    1.0        4 mm
#> 5  1.02748190         FALSE           0           3.375    1.0        8 mm
#> 6 -0.06587817         FALSE           0           3.375    1.0        8 mm
#> 7  0.00000000         FALSE           0        1728.000    1.0       voxel
#> 8  4.09490015         FALSE           0        1728.000    0.8       voxel
#>                                                 curve_use
#> 1                               spatial-frequency readout
#> 2 localized squared evidence; not a spectral-energy curve
#> 3                               spatial-frequency readout
#> 4 localized squared evidence; not a spectral-energy curve
#> 5                               spatial-frequency readout
#> 6 localized squared evidence; not a spectral-energy curve
#> 7                               spatial-frequency readout
#> 8 localized squared evidence; not a spectral-energy curve
```

## What if you only have a statistic map?

You do not need to store all permutation maps in memory. If an upstream
workflow can regenerate the `b`th null statistic map, pass that
deterministic draw function directly.

``` r

draw_null_map <- function(b) {
  glm_null$draw(b)
}

fit_from_function <- scale_space_energy(
  stat_vol = glm_null$stat_vol,
  mask = glm_null$mask,
  null_fun = draw_null_map,
  n_perm = glm_null$n_perm,
  scales = scale_space(c(4, 8), include_global = TRUE, include_identity = TRUE),
  statistics = list(global_lowpass_energy(), global_band_energy(), local_energy())
)

scale_space_diagnostics(fit_from_function)[, c(
  "channel", "frequency_readout", "best_scale", "best_p_fwer", "curve_use"
)]
#>                 channel frequency_readout best_scale best_p_fwer
#> 1    global_band_energy              TRUE     global         1.0
#> 2 global_lowpass_energy              TRUE     global         1.0
#> 3          local_energy             FALSE      voxel         0.8
#>                                                 curve_use
#> 1                               spatial-frequency readout
#> 2                               spatial-frequency readout
#> 3 localized squared evidence; not a spectral-energy curve
```

The draw function must be deterministic by index: `draw_null_map(7)`
should return the same map every time.
[`scale_space_energy()`](https://bbuchsbaum.github.io/neurothresh/reference/scale_space_energy.md)
makes multiple passes over the same null family to estimate null moments
and then compute the max-null distribution.

``` r

external_null <- null_from_function(
  draw = draw_null_map,
  stat_vol = glm_null$stat_vol,
  mask = glm_null$mask,
  n_perm = glm_null$n_perm,
  kind = "external_stat_map_generator"
)

fit_from_null_spec <- scale_space_energy(
  external_null,
  scales = scale_space(c(4, 8), include_global = TRUE, include_identity = TRUE),
  statistics = list(global_band_energy(), local_energy())
)

summary(fit_from_null_spec)[, c("scale_label", "channel", "p_fwer", "curve_use")]
#>   scale_label            channel p_fwer
#> 1      global global_band_energy    1.0
#> 2      global       local_energy    1.0
#> 3        4 mm global_band_energy    1.0
#> 4        4 mm       local_energy    1.0
#> 5        8 mm global_band_energy    1.0
#> 6        8 mm       local_energy    1.0
#> 7       voxel global_band_energy    1.0
#> 8       voxel       local_energy    0.8
#>                                                 curve_use
#> 1                               spatial-frequency readout
#> 2 localized squared evidence; not a spectral-energy curve
#> 3                               spatial-frequency readout
#> 4 localized squared evidence; not a spectral-energy curve
#> 5                               spatial-frequency readout
#> 6 localized squared evidence; not a spectral-energy curve
#> 7                               spatial-frequency readout
#> 8 localized squared evidence; not a spectral-energy curve
```

If you already have precomputed null maps, use
[`null_from_maps()`](https://bbuchsbaum.github.io/neurothresh/reference/null_from_maps.md)
instead. The statistical requirement is the same in both cases: the
observed map and every null map must be produced by the same upstream
statistic and must pass through the same scale-space transforms before
correction.

``` r

stored_null_maps <- lapply(seq_len(glm_null$n_perm), draw_null_map)
maps_null <- null_from_maps(
  stat_vol = glm_null$stat_vol,
  null_maps = stored_null_maps,
  mask = glm_null$mask
)

fit_from_maps <- scale_space_energy(
  maps_null,
  scales = scale_space(8, include_global = TRUE, include_identity = TRUE),
  statistics = list(global_band_energy(), local_energy())
)

scale_space_diagnostics(fit_from_maps)[, c("channel", "best_scale", "best_p_fwer")]
#>              channel best_scale best_p_fwer
#> 1 global_band_energy     global         1.0
#> 2       local_energy      voxel         0.8
```

This approach is related to post-statistic spatial filtering methods
such as LISA: the same spatial operation is applied to observed and
permuted maps before inference. The distinction is that scale-space
localized energy inference searches a prespecified resolution family and
uses max-null FWER over the full search. It is therefore not smoothing
until significance; the smoothing scales are part of the corrected test
family.

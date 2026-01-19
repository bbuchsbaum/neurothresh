# Estimate Gaussian ACF FWHM from a single map (map-only)

Map-only smoothness estimator based on robust sign correlations at short
lags. This assumes approximate stationarity and that signal occupies a
relatively small fraction of the mask.

## Usage

``` r
estimate_fwhm_vox(
  z_vol,
  mask = NULL,
  clip = 3,
  isotropic = FALSE,
  z_trim = 1.5
)
```

## Arguments

- z_vol:

  3D volume-like object.

- mask:

  Optional logical mask aligned to `z_vol`. If NULL, uses
  `is.finite(as.numeric(z_vol))`.

- clip:

  Clip value applied to `z_vol` prior to taking `sign`, to reduce
  sensitivity to large outliers (default 3).

- isotropic:

  Logical; if TRUE, returns a single averaged FWHM replicated over
  x/y/z.

- z_trim:

  If not NULL, restricts estimation to voxels with `|Z| <= z_trim`
  (default 1.5) as a heuristic to reduce contamination from true signal
  when only a single statistic map is available.

## Value

Numeric vector length 3: estimated `fwhm_vox` in voxel units.

## Details

This estimates an effective stationary smoothness parameter from the
statistic field itself (not from residuals). It is intended for
model-based Monte Carlo calibration when subject-level/residual
resampling is not available.

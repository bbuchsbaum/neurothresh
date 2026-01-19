# RFT Peak-Level FWER

Computes the FWER-corrected threshold for peak (voxelwise) inference
using Random Field Theory.

## Usage

``` r
rft_peak_fwer(
  stat_vol,
  mask = NULL,
  fwhm_mm,
  alpha = 0.05,
  df = Inf,
  tail = c("pos", "neg", "two")
)
```

## Arguments

- stat_vol:

  NeuroVol containing statistic map

- mask:

  Optional LogicalNeuroVol for analysis domain

- fwhm_mm:

  Smoothness in mm. Scalar or length-3 vector for anisotropic
  smoothness.

- alpha:

  Significance level (default 0.05)

- df:

  Degrees of freedom. Use Inf for Gaussian field (Z-scores).

- tail:

  "pos" (activations), "neg" (deactivations), or "two"

## Value

List with components:

- threshold:

  The FWER-corrected threshold

- sig_mask:

  Logical mask of significant voxels

- method:

  "rft_peak"

## Details

Uses the expected Euler characteristic approximation: \$\$P(\max Z \> u)
\approx \sum\_{d=0}^{D} \text{Resels}\_d \cdot EC_d(u)\$\$

For high thresholds in smooth Gaussian fields, this provides accurate
FWER control.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- rft_peak_fwer(z_map, fwhm_mm = 8, alpha = 0.05)
} # }
```

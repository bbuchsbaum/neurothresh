# RFT Peak-Level FWER

Computes the FWER-corrected threshold for peak (voxelwise) inference
using Random Field Theory.

## Usage

``` r
rft_peak_fwer(
  stat_vol,
  mask = NULL,
  alpha = 0.05,
  tail = c("pos", "neg", "two"),
  smoothness,
  df = Inf
)
```

## Arguments

- stat_vol:

  3D Z-equivalent statistic map (array or `NeuroVol`). Pass
  already-canonicalized input (use
  [`canonicalize()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize.md)
  first for t-statistics or -log10(p) maps).

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- alpha:

  Significance level (default 0.05).

- tail:

  One of `"pos"`, `"neg"`, `"two"`. For `"two"`, each direction uses
  `alpha / 2` and the union is reported, targeting FWER `alpha` across
  both tails.

- smoothness:

  A
  [`smoothness_fwhm()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_fwhm.md)
  object specifying field smoothness in millimetres. Numeric input is
  also accepted (length 1 or 3, mm) and coerced internally.

- df:

  Degrees of freedom. Use `Inf` for Gaussian field (Z-scores).

## Value

A
[`neurothresh_result`](https://bbuchsbaum.github.io/neurothresh/reference/new_neurothresh_result.md)
with `method = "rft_peak_fwer"`. Slots: `threshold` is the positive-side
critical value; `reject` is the 3D logical rejection volume; `pvals` is
`NA` (RFT does not produce per-voxel corrected p-values without further
computation).

## Details

Uses the expected Euler characteristic approximation: \$\$P(\max Z \> u)
\approx \sum\_{d=0}^{D} \mathrm{Resels}\_d \cdot EC_d(u)\$\$

For high thresholds in smooth Gaussian fields, this provides accurate
FWER control.

## See also

[`rft_cluster_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_cluster_fwer.md)
for cluster-level inference;
[`tfce_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md)
and
[`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md)
for alternative methods.

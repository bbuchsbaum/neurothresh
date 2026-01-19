# TFCE Transform

Computes the Threshold-Free Cluster Enhancement transform: \$\$TFCE(v) =
\int_0^{h_v} e(h)^E \cdot h^H \\ dh\$\$

## Usage

``` r
tfce_transform(
  stat_vol,
  mask = NULL,
  H = 2,
  E = 0.5,
  dh = 0.1,
  tail = c("pos", "neg")
)
```

## Arguments

- stat_vol:

  NeuroVol containing statistic map

- mask:

  Optional LogicalNeuroVol for analysis domain

- H:

  Height exponent (default 2.0)

- E:

  Extent exponent (default 0.5)

- dh:

  Threshold increment (default 0.1)

- tail:

  "pos" for activations, "neg" for deactivations

## Value

NeuroVol with TFCE-transformed values

## Details

where e(h) is the extent of the cluster containing voxel v at threshold
h.

TFCE provides threshold-free inference by integrating cluster extent
over a range of thresholds. Default parameters (H=2, E=0.5) are those
recommended by Smith & Nichols (2009).

## References

Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
enhancement: Addressing problems of smoothing, threshold dependence and
localisation in cluster inference. NeuroImage, 44(1), 83-98.

## Examples

``` r
if (FALSE) { # \dontrun{
tfce_map <- tfce_transform(z_map, H = 2, E = 0.5)
} # }
```

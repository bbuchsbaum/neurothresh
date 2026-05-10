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
  tail = c("pos", "neg", "two")
)
```

## Arguments

- stat_vol:

  3D Z-equivalent statistic map (array or `NeuroVol`).

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- H:

  Height exponent (default 2.0).

- E:

  Extent exponent (default 0.5).

- dh:

  Threshold increment (default 0.1).

- tail:

  One of `"pos"`, `"neg"`, `"two"`. For `"two"` the result is the
  per-voxel maximum of the positive- and negative-tail TFCE maps.

## Value

Volume (same type as `stat_vol`) with TFCE-transformed values.

## Details

where \\e(h)\\ is the extent of the cluster containing voxel \\v\\ at
threshold \\h\\.

## References

Smith & Nichols (2009), NeuroImage 44(1), 83-98.

## See also

[`tfce_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md).

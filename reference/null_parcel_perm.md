# Parcel-level permutation null

Used internally by
[`hier_scan()`](https://bbuchsbaum.github.io/neurothresh/reference/hier_scan.md).
For each draw, voxel signs are flipped independently at the voxel level
(the same family used by
[`null_voxel_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_voxel_signflip.md)
but exposed for hierarchical workflows).

## Usage

``` r
null_parcel_perm(
  stat_vol,
  parcels = NULL,
  mask = NULL,
  n_perm = 5000,
  seed = NULL
)
```

## Arguments

- stat_vol:

  Observed (canonicalized) Z-equivalent map.

- parcels:

  A `ClusteredNeuroVol`, integer label array, or vector aligned with
  `stat_vol`.

- mask:

  Optional logical 3D mask.

- n_perm:

  Number of draws.

- seed:

  Optional random seed.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
with `kind = "parcel_perm"`.

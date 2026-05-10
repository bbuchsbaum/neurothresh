# Voxel-level sign-flip null

Permutes the sign of each voxel independently. This null is
anti-conservative in the presence of spatial correlation; prefer
subject-level sign-flips
([`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md)
/
[`null_subject_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_subject_signflip.md))
when subject data are available.

## Usage

``` r
null_voxel_signflip(stat_vol, mask = NULL, n_perm = 5000, seed = NULL)
```

## Arguments

- stat_vol:

  Observed (canonicalized) Z-equivalent statistic map.

- mask:

  Optional logical 3D mask. Defaults to finite voxels in `stat_vol`.

- n_perm:

  Number of sign-flip draws.

- seed:

  Optional random seed; the resulting `draw(b)` is deterministic by
  index.

## Value

A
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md).

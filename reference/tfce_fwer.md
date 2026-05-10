# TFCE with Permutation FWER

Computes TFCE with permutation-based FWER control using the maxT
procedure. The null distribution is supplied via a
[`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md);
for subject-level inference build one with
[`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md)
/
[`null_two_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_two_sample.md)
/
[`null_glm()`](https://bbuchsbaum.github.io/neurothresh/reference/null_glm.md).
For a fast voxel-level fallback use
[`null_voxel_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_voxel_signflip.md)
(note: voxel sign-flips are calibrated only when the field is
uncorrelated under the null).

## Usage

``` r
tfce_fwer(
  stat_vol,
  mask = NULL,
  alpha = 0.05,
  tail = c("pos", "neg", "two"),
  H = 2,
  E = 0.5,
  dh = 0.1,
  null_spec = NULL,
  null_fun = NULL,
  n_perm = NULL,
  allow_voxel_signflip = FALSE,
  seed = NULL
)
```

## Arguments

- stat_vol:

  3D Z-equivalent statistic map.

- mask:

  Optional logical 3D mask.

- alpha:

  Significance level (default 0.05).

- tail:

  One of `"pos"`, `"neg"`, `"two"`.

- H:

  Height exponent (default 2.0).

- E:

  Extent exponent (default 0.5).

- dh:

  Threshold increment (default 0.1).

- null_spec:

  A
  [`null_spec`](https://bbuchsbaum.github.io/neurothresh/reference/null_spec.md)
  supplying `n_perm` deterministic null draws.

- null_fun:

  Optional function `function(b)` returning a null draw for permutation
  `b`. Kept for compatibility with the package's older permutation API;
  prefer `null_spec` for new code.

- n_perm:

  Number of permutations to consume from `null_spec`. Defaults to
  `null_spec$n_perm`.

- allow_voxel_signflip:

  Logical; if `TRUE`, allows a fallback to independent voxel sign-flips
  when neither `null_spec` nor `null_fun` is supplied. This fallback is
  fast but can be anti-conservative for spatially correlated fields.

- seed:

  Optional random seed; affects only RNG state outside `null_spec` (the
  spec itself is deterministic by index).

## Value

A
[`neurothresh_result`](https://bbuchsbaum.github.io/neurothresh/reference/new_neurothresh_result.md)
with `method = "tfce_fwer"`. `pvals` is the voxelwise FWER-corrected
p-value volume; `extras$tfce` holds the observed TFCE map; `null_max` is
the vector of per-permutation max-TFCE values.

## See also

[`tfce_transform()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_transform.md),
[`null_voxel_signflip()`](https://bbuchsbaum.github.io/neurothresh/reference/null_voxel_signflip.md),
[`null_one_sample()`](https://bbuchsbaum.github.io/neurothresh/reference/null_one_sample.md).

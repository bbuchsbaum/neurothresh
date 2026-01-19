# Alpha-spending hierarchical descent without step-down

Tests each node against its own permutation distribution and only
descends when the node is significant at its allocated alpha.

## Usage

``` r
hier_descend_alpha(
  idx_parent,
  bbox,
  alpha_budget,
  z_vec,
  pi_vec,
  x,
  y,
  z,
  perm_fun,
  n_perm,
  kappa_grid,
  gamma = 0.5,
  min_alpha = 1e-06,
  min_voxels = 8
)
```

## Arguments

- idx_parent:

  Integer vector of parent mask-space indices

- bbox:

  Integer vector of length 6 for parent bounding box

- alpha_budget:

  Alpha budget for this subtree

- z_vec:

  Numeric vector of Z-scores (mask-space)

- pi_vec:

  Numeric vector of prior weights (mask-space)

- x, y, z:

  Integer vectors of coordinates (mask-space)

- perm_fun:

  Function that takes permutation index and returns Z vector

- n_perm:

  Number of permutations

- kappa_grid:

  Numeric vector of positive kappa values

- gamma:

  Fraction of alpha for testing at each level (default 0.5)

- min_alpha:

  Minimum alpha to continue testing (default 1e-6)

- min_voxels:

  Minimum region size to test (default 8)

## Value

List of discovered regions with scores and p-values

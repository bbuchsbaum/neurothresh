# Score all children in one pass through parent voxels

Computes the omnibus score (max of U_0 and best S_kappa) for all
children in a single O(n) pass through the parent voxels.

## Usage

``` r
score_children_onepass(Z_vec, idx_parent, split, kappa_grid, do_abs = FALSE)
```

## Arguments

- Z_vec:

  Numeric vector of Z-scores (mask-space, length N)

- idx_parent:

  Integer vector of parent mask-space indices

- split:

  List returned by
  [`octree_split_info`](https://bbuchsbaum.github.io/neurothresh/reference/octree_split_info.md)

- kappa_grid:

  Numeric vector of positive kappa values

- do_abs:

  Logical, whether to take absolute value of Z (for two-sided)

## Value

Numeric vector of length m with child scores

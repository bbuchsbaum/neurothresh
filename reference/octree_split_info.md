# Build octree split info for a parent node

Computes child assignments and caches constants needed for one-pass
scoring of all children.

## Usage

``` r
octree_split_info(idx_parent, bbox, x, y, z, pi_vec, min_pi_mass = 1e-10)
```

## Arguments

- idx_parent:

  Integer vector of 1-based mask-space indices

- bbox:

  Integer vector of length 6: (x0, x1, y0, y1, z0, z1)

- x:

  Integer vector of x coordinates (mask-space, length N)

- y:

  Integer vector of y coordinates (mask-space, length N)

- z:

  Integer vector of z coordinates (mask-space, length N)

- pi_vec:

  Numeric vector of prior weights (mask-space, length N)

- min_pi_mass:

  Minimum prior mass to keep a child (default 1e-10)

## Value

A list with components:

- m:

  Number of non-empty children (0 to 8)

- child_id:

  Integer vector of child assignments (-1 or 0..m-1)

- w_parent:

  Numeric vector of prior weights for idx_parent

- log_den1:

  Numeric vector of log(sum(pi)) per child

- inv_sqrt_den2:

  Numeric vector of 1/sqrt(sum(pi^2)) per child

- keep8:

  Integer vector indicating which of 1..8 octants are kept

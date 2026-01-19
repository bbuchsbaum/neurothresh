# Split indices with bounding boxes and prior masses

Like
[`split_indices_bbox`](https://bbuchsbaum.github.io/neurothresh/reference/split_indices_bbox.md)
but also computes sum(pi) and sum(pi^2) for each child for alpha
allocation.

## Usage

``` r
split_indices_bbox_mass(idx_parent, split, x, y, z, pi_vec)
```

## Arguments

- idx_parent:

  Integer vector of parent mask-space indices

- split:

  List returned by
  [`octree_split_info`](https://bbuchsbaum.github.io/neurothresh/reference/octree_split_info.md)

- x:

  Integer vector of x coordinates

- y:

  Integer vector of y coordinates

- z:

  Integer vector of z coordinates

- pi_vec:

  Numeric vector of prior weights

## Value

List with components:

- idx:

  List of m IntegerVectors with child indices

- bbox:

  Integer matrix (m x 6) with child bounding boxes

- den1:

  Numeric vector of sum(pi) per child

- den2:

  Numeric vector of sum(pi^2) per child

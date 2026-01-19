# Split indices with bounding boxes

Splits parent indices into children and computes bounding box for each
child in a single pass.

## Usage

``` r
split_indices_bbox(idx_parent, split, x, y, z)
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

## Value

List with components:

- idx:

  List of m IntegerVectors with child indices

- bbox:

  Integer matrix (m x 6) with child bounding boxes

# Split parent indices into child index lists

Partitions the parent voxel indices into lists for each child, for
recursive descent.

## Usage

``` r
split_child_indices(idx_parent, split)
```

## Arguments

- idx_parent:

  Integer vector of parent mask-space indices

- split:

  List returned by
  [`octree_split_info`](https://bbuchsbaum.github.io/neurothresh/reference/octree_split_info.md)

## Value

List of m IntegerVectors, one per child

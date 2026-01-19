# Recursive octree split

Recursively partitions a region into octree children until regions are
smaller than min_voxels.

## Usage

``` r
octree_split(indices, x, y, z, min_voxels = 8)
```

## Arguments

- indices:

  Integer vector of mask-space indices

- x:

  Integer vector of x coordinates (mask-space)

- y:

  Integer vector of y coordinates (mask-space)

- z:

  Integer vector of z coordinates (mask-space)

- min_voxels:

  Minimum region size to continue splitting

## Value

List of index vectors for leaf regions

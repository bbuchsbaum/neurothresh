# Compute bounding box from mask-space indices

Compute bounding box from mask-space indices

## Usage

``` r
bbox_from_indices(idx, x, y, z)
```

## Arguments

- idx:

  Integer vector of 1-based mask-space indices

- x:

  Integer vector of x coordinates (mask-space)

- y:

  Integer vector of y coordinates (mask-space)

- z:

  Integer vector of z coordinates (mask-space)

## Value

Integer vector of length 6: (x0, x1, y0, y1, z0, z1)

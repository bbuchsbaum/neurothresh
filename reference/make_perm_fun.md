# Create permutation function for hierarchical analysis

Creates a function that returns sign-flipped Z vectors for use in
hierarchical descent.

## Usage

``` r
make_perm_fun(z_vec, n_perm, seed = NULL)
```

## Arguments

- z_vec:

  Numeric vector of Z-scores (mask-space)

- n_perm:

  Number of permutations

- seed:

  Random seed for reproducibility (optional)

## Value

A function that takes permutation index b (1 to n_perm) and returns a
sign-flipped Z vector

## Details

The sign matrix is pre-generated and stored, so calling the returned
function multiple times with the same index gives the same result.

# Generate null scores for parcel-level analysis

Vectorized null score generation for parcels using efficient group-wise
operations.

## Usage

``` r
generate_null_scores_parcels(
  z_vec,
  pi_vec,
  lab_vec,
  n_perm,
  kappa_grid = c(0.5, 1, 2),
  seed = NULL
)
```

## Arguments

- z_vec:

  Numeric vector of Z-scores (mask-space)

- pi_vec:

  Numeric vector of prior weights (mask-space)

- lab_vec:

  Integer vector of parcel labels (mask-space)

- n_perm:

  Number of permutations

- kappa_grid:

  Numeric vector of positive kappa values

- seed:

  Random seed for reproducibility (optional)

## Value

Matrix of dimension (n_perm x n_parcels) with null scores

## Details

This is more efficient than
[`generate_null_scores`](https://bbuchsbaum.github.io/neurothresh/reference/generate_null_scores.md)
for parcel-level analysis because it uses vectorized group operations
(rowsum) instead of looping over regions.

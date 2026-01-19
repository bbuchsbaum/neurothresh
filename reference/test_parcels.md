# Test at parcel level with step-down

Tests all parcels simultaneously using Westfall-Young step-down, then
returns significant parcels for further descent.

## Usage

``` r
test_parcels(
  z_vec,
  pi_vec,
  lab_vec,
  n_perm,
  kappa_grid,
  alpha = 0.05,
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

- alpha:

  Significance level

- seed:

  Random seed (optional)

## Value

List with components:

- stepdown:

  Data frame from wy_stepdown

- sig_parcels:

  Integer vector of significant parcel labels

- parcel_indices:

  List of index vectors for each parcel

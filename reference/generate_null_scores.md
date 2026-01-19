# Generate null distribution via sign-flipping

Generates null scores for regions by randomly sign-flipping the Z-score
map. This is valid when the original data are symmetric around zero
under the null (e.g., group difference maps).

## Usage

``` r
generate_null_scores(
  z_vec,
  pi_vec,
  regions,
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

- regions:

  List of index vectors defining regions

- n_perm:

  Number of permutations

- kappa_grid:

  Numeric vector of positive kappa values

- seed:

  Random seed for reproducibility (optional)

## Value

Matrix of dimension (n_perm x n_regions) with null scores

## Details

For each permutation, a random sign (+1 or -1) is drawn independently
for each voxel and multiplied with the Z-scores. The omnibus score (max
of U_0 and best S_kappa) is computed for each region.

## Examples

``` r
if (FALSE) { # \dontrun{
z_vec <- rnorm(1000)
pi_vec <- rep(1, 1000)
regions <- list(1:100, 101:200, 201:300)
null_matrix <- generate_null_scores(z_vec, pi_vec, regions,
                                    n_perm = 1000)
} # }
```

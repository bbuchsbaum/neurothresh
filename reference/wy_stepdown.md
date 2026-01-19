# Westfall-Young Step-Down Procedure

Implements the Westfall-Young step-down maxT procedure for strong
control of the family-wise error rate (FWER).

## Usage

``` r
wy_stepdown(observed, null_matrix, alpha = 0.05)
```

## Arguments

- observed:

  Numeric vector of observed test statistics

- null_matrix:

  Matrix of null statistics (n_perm x n_tests)

- alpha:

  Significance level (default 0.05)

## Value

A data.frame with columns:

- region:

  Index of the region (1 to n_tests)

- score:

  Observed test statistic

- p_adj:

  Adjusted p-value (monotonically non-decreasing)

- rejected:

  Logical, whether the null is rejected

## Details

The Westfall-Young step-down procedure provides strong FWER control with
higher power than single-step methods by:

1.  Sorting regions by observed score (descending)

2.  Computing successive maxima from the null distribution

3.  Ensuring adjusted p-values are monotonically non-decreasing

This is the recommended method for hierarchical LR-MFT as it preserves
the logical ordering of rejections.

## References

Westfall, P. H., & Young, S. S. (1993). Resampling-based multiple
testing: Examples and methods for p-value adjustment. Wiley.

## Examples

``` r
if (FALSE) { # \dontrun{
# Observed scores for 5 regions
observed <- c(3.5, 2.1, 4.2, 1.8, 2.9)

# Null distribution from 1000 permutations
null_matrix <- matrix(rnorm(5000), nrow = 1000, ncol = 5)

# Run step-down procedure
result <- wy_stepdown(observed, null_matrix, alpha = 0.05)
} # }
```

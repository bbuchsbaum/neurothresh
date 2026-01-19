# Single-step MaxT Procedure

Implements the single-step maxT procedure for FWER control. Less
powerful than step-down but simpler.

## Usage

``` r
maxT_singlestep(observed, null_matrix, alpha = 0.05)
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

  Index of the region

- score:

  Observed test statistic

- p_adj:

  Adjusted p-value

- rejected:

  Logical, whether the null is rejected

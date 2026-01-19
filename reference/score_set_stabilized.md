# Compute variance-stabilized score U_0(R)

Computes the variance-stabilized diffuse detection score: \$\$U_0(R) =
\frac{\sum\_{v \in R} \pi(v) \cdot Z(v)}{\sqrt{\sum\_{v \in R}
\pi(v)^2}}\$\$

## Usage

``` r
score_set_stabilized(indices, z_vec, pi_vec)
```

## Arguments

- indices:

  Integer vector of 1-based mask-space indices defining the region

- z_vec:

  Numeric vector of Z-scores in mask-space

- pi_vec:

  Numeric vector of prior weights in mask-space

## Value

A list with components:

- U0:

  The variance-stabilized score

- n_eff:

  The effective sample size: (sum(pi))^2 / sum(pi^2)

## Details

This provides fair comparison across regions of different sizes by
normalizing by the effective sample size.

Under the null hypothesis, U_0(R) has approximately unit variance
regardless of region size, enabling fair comparison across scales.

The effective sample size `n_eff` represents the equivalent number of
independent observations.

## See also

[`score_set`](https://bbuchsbaum.github.io/neurothresh/reference/score_set.md)
for kappa \> 0 focal detection

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare scores across regions of different sizes
z_vec <- rnorm(1000)
pi_vec <- rep(1, 1000)

# Small region
small <- score_set_stabilized(1:10, z_vec, pi_vec)

# Large region
large <- score_set_stabilized(1:100, z_vec, pi_vec)

# Both have comparable variance under null
} # }
```

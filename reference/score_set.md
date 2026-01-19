# Compute prior-weighted soft evidence statistic T_kappa(R)

Computes the prior-weighted log-sum-exp statistic for a region R:
\$\$T\_\kappa(R) = \log \sum\_{v \in R} \pi(v) \cdot \exp(\kappa \cdot
Z(v))\$\$

## Usage

``` r
score_set(indices, z_vec, pi_vec, kappa = 1)
```

## Arguments

- indices:

  Integer vector of 1-based mask-space indices defining the region

- z_vec:

  Numeric vector of Z-scores in mask-space

- pi_vec:

  Numeric vector of prior weights in mask-space

- kappa:

  Numeric scalar, temperature parameter (default 1.0). Higher kappa
  emphasizes peaks; kappa -\> 0 gives diffuse detection.

## Value

Numeric scalar, the T_kappa score for the region

## Details

This is the core statistic for LR-MFT. For kappa = 0, use
[`score_set_stabilized`](https://bbuchsbaum.github.io/neurothresh/reference/score_set_stabilized.md)
instead for variance-stabilized diffuse detection.

## See also

[`score_set_stabilized`](https://bbuchsbaum.github.io/neurothresh/reference/score_set_stabilized.md)
for kappa = 0 case

## Examples

``` r
if (FALSE) { # \dontrun{
# Score a region with 100 voxels
z_vec <- rnorm(1000)
pi_vec <- rep(1, 1000)
indices <- 1:100
score <- score_set(indices, z_vec, pi_vec, kappa = 1.0)
} # }
```

# Compute omnibus score combining diffuse and focal detection

Computes the maximum of the variance-stabilized diffuse score U_0 and
the best soft-max score S_kappa across a grid of kappa values.

## Usage

``` r
score_set_omnibus(indices, z_vec, pi_vec, kappa_grid = c(0.5, 1, 2))
```

## Arguments

- indices:

  Integer vector of 1-based mask-space indices

- z_vec:

  Numeric vector of Z-scores in mask-space

- pi_vec:

  Numeric vector of prior weights in mask-space

- kappa_grid:

  Numeric vector of positive kappa values to try

## Value

A list with components:

- T_omnibus:

  The maximum score across all statistics

- U0:

  The variance-stabilized diffuse score

- S_best:

  The best soft-max score

- kappa_best:

  The kappa value that gave S_best

## Details

The omnibus statistic is: \$\$T(R) = \max(U_0(R), \max\_{\kappa}
S\_\kappa(R))\$\$

where S_kappa is the normalized soft-max: \$\$S\_\kappa(R) =
\frac{1}{\kappa} \left\[ \log \sum \pi \exp(\kappa Z) - \log \sum \pi
\right\]\$\$

# AFNI-style ACF smoothness specification

Tagged container for the three-parameter mixed Gaussian +
monoexponential autocorrelation model used by AFNI's `3dFWHMx -acf` and
downstream Monte-Carlo null generators.

## Usage

``` r
smoothness_acf(a, b, c)
```

## Arguments

- a:

  Numeric scalar in (0, 1); mixing weight of the Gaussian component.

- b:

  Numeric scalar \> 0; Gaussian width parameter (mm).

- c:

  Numeric scalar \> 0; exponential length scale (mm).

## Value

Object of class
`c("neurothresh_smoothness_acf", "neurothresh_smoothness")`.

## Examples

``` r
s <- smoothness_acf(a = 0.7, b = 3.0, c = 8.0)
```

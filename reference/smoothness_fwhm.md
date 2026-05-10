# FWHM-based smoothness specification

Tagged container for full-width-at-half-maximum (FWHM) smoothness in
millimetres. RFT methods accept only this form; Monte-Carlo null
factories may also accept it (alongside
[`smoothness_acf()`](https://bbuchsbaum.github.io/neurothresh/reference/smoothness_acf.md)).

## Usage

``` r
smoothness_fwhm(fwhm_mm)
```

## Arguments

- fwhm_mm:

  Numeric; either a length-1 isotropic FWHM or a length-3 per-axis FWHM,
  in millimetres. Must be strictly positive.

## Value

Object of class
`c("neurothresh_smoothness_fwhm", "neurothresh_smoothness")`.

## Examples

``` r
s1 <- smoothness_fwhm(6)
s3 <- smoothness_fwhm(c(6, 6, 8))
```

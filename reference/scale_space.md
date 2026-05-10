# Define a Gaussian Scale-Space Grid

Define a Gaussian Scale-Space Grid

## Usage

``` r
scale_space(
  fwhm_mm,
  include_identity = TRUE,
  include_global = TRUE,
  report = c("fwhm", "effective_basis")
)
```

## Arguments

- fwhm_mm:

  Numeric vector of Gaussian FWHM scales in millimeters.

- include_identity:

  Include the unsmoothed voxelwise endpoint.

- include_global:

  Include the global endpoint.

- report:

  Metadata requested by summaries.

## Value

A `scale_space_spec` object.

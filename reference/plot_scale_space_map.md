# Extract a Scale-Space Corrected Map

Extract a Scale-Space Corrected Map

## Usage

``` r
plot_scale_space_map(
  fit,
  scale_index = NULL,
  channel = NULL,
  what = c("p_corr", "statistic", "raw")
)
```

## Arguments

- fit:

  A `scale_space_energy_result`.

- scale_index:

  Scale index to extract.

- channel:

  Channel name to extract.

- what:

  One of "p_corr", "statistic", or "raw".

## Value

A volume-like object aligned to the input statistic map.

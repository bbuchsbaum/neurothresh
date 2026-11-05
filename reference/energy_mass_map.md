# Above-Null Energy Mass Map

Above-Null Energy Mass Map

## Usage

``` r
energy_mass_map(fit, scale_index = NULL, mass = NULL)
```

## Arguments

- fit:

  A `scale_space_energy_result`.

- scale_index:

  Scale index to extract. Defaults to the strongest finite
  `local_energy` scale. The global endpoint is not valid for an energy
  mass map because it is spatially constant.

- mass:

  Optional cumulative mass cutoff.

## Value

A volume-like normalized above-null energy mass map.

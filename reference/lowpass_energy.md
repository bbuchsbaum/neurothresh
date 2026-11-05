# Low-Pass Energy Statistic Channel

Computes energy after smoothing the signed statistic field. With
`scope = "global"`, the statistic is \\\sum_v \[S\_\rho z(v)\]^2\\, a
cumulative low-frequency energy curve. With `scope = "local"`, the
statistic is \\\[S\_\rho z(v)\]^2\\.

## Usage

``` r
lowpass_energy(scope = c("global", "local"))
```

## Arguments

- scope:

  `"global"` returns one whole-map statistic per scale; `"local"`
  returns a location-wise field.

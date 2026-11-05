# Band-Pass Energy Statistic Channel

Explicit opt-in spatial band-pass channel using a difference of
smoothers around the requested scale. The local form returns a field
\\\[S\_{\rho/r}z(v) - S\_{\rho r}z(v)\]^2\\. The global form returns
\\\sum_v \[S\_{\rho/r}z(v) - S\_{\rho r}z(v)\]^2\\, which is the
clearest curve for reading out the dominant spatial-frequency band.

## Usage

``` r
bandpass_energy(ratio = sqrt(2), scope = c("local", "global"))
```

## Arguments

- ratio:

  Multiplicative scale ratio around the target FWHM.

- scope:

  `"local"` preserves the historical location-wise behavior; `"global"`
  returns one whole-map statistic per scale.

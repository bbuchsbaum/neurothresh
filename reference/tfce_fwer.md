# TFCE with Permutation FWER

Computes TFCE with permutation-based FWER control using the maxT
procedure.

## Usage

``` r
tfce_fwer(
  stat_vol,
  mask = NULL,
  n_perm = 5000,
  H = 2,
  E = 0.5,
  dh = 0.1,
  alpha = 0.05,
  seed = NULL,
  tail = c("pos", "neg", "two")
)
```

## Arguments

- stat_vol:

  NeuroVol containing statistic map

- mask:

  Optional LogicalNeuroVol for analysis domain

- n_perm:

  Number of permutations (default 5000)

- H:

  Height exponent (default 2.0)

- E:

  Extent exponent (default 0.5)

- dh:

  Threshold increment (default 0.1)

- alpha:

  Significance level (default 0.05)

- seed:

  Random seed for reproducibility

- tail:

  "pos", "neg", or "two"

## Value

List with components:

- tfce:

  TFCE-transformed observed map

- threshold:

  FWER-corrected threshold

- sig_mask:

  Logical mask of significant voxels

- p_corr:

  Voxelwise FWER-corrected p-values

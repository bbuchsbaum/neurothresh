# neurothresh: Neuroimaging Statistical Thresholding with LR-MFT

Statistical thresholding for neuroimaging data using Likelihood-Ratio
Matched-Filter Thresholding (LR-MFT) with hierarchical stepdown
inference. Provides theoretically superior alternative to traditional
GRF/RFT methods for detecting spatially extended activations. Includes
baseline implementations (RFT, TFCE, cluster-FDR) for comparison.

Statistical thresholding for neuroimaging data using Likelihood-Ratio
Matched-Filter Thresholding (LR-MFT) with hierarchical stepdown
inference.

## Main Functions

- [`hier_scan`](https://bbuchsbaum.github.io/neurothresh/reference/hier_scan.md):

  Main entry point for hierarchical LR-MFT analysis

- [`score_set`](https://bbuchsbaum.github.io/neurothresh/reference/score_set.md):

  Compute prior-weighted statistic T_kappa(R)

- [`canonicalize_stat`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize_stat.md):

  Convert Z/t/neglog10p to Z-scores

## Baseline Methods

- [`rft_peak_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/rft_peak_fwer.md):

  RFT peak-level FWER

- [`rft_cluster_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/rft_cluster_fwer.md):

  RFT cluster-level FWER

- [`tfce_transform`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_transform.md):

  TFCE transformation

- [`tfce_fwer`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md):

  TFCE with permutation FWER

- [`cluster_fdr`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md):

  Topological cluster-FDR

## Key Concepts

LR-MFT is provably more powerful than traditional GRF/RFT methods for
detecting spatially extended activation blobs. It uses a prior-weighted
soft evidence statistic:

\$\$T\_\kappa(R) = \log \sum\_{v \in R} \pi(v) \cdot \exp(\kappa \cdot
Z(v))\$\$

where \\\pi(v)\\ is the prior weight and \\\kappa\\ is the temperature.

## See also

Useful links:

- <https://github.com/bbuchsbaum/neurothresh>

- <https://bbuchsbaum.github.io/neurothresh/>

- Report bugs at <https://github.com/bbuchsbaum/neurothresh/issues>

## Author

**Maintainer**: Brad Buchsbaum <brad.buchsbaum@gmail.com>

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

**neurothresh** is an R package for neuroimaging statistical
thresholding that implements: 1. **LR-MFT (Likelihood-Ratio
Matched-Filter Thresholding)** with hierarchical stepdown inference - a
prior-weighted soft evidence aggregation method that is provably more
powerful than RFT peak-height for multi-focal/distributed brain
activations 2. **Standard baseline methods**: RFT peak/cluster FWER,
TFCE, cluster-FDR

The package integrates with **neuroim2** (volumetric brain imaging data
structures) and **neuroatlas** (brain parcellations).

## Build & Development Commands

``` bash
# Install dependencies
Rscript -e "install.packages(c('Rcpp', 'testthat', 'roxygen2', 'devtools'))"
Rscript -e "devtools::install_deps()"

# Build package
Rscript -e "devtools::document()"
R CMD INSTALL .

# Run all tests
Rscript -e "devtools::test()"

# Run single test file
Rscript -e "testthat::test_file('tests/testthat/test-score-node.R')"

# Check package (CRAN-style)
R CMD check --as-cran .

# Load for interactive development
Rscript -e "devtools::load_all()"

# Compile Rcpp attributes after modifying C++ code
Rscript -e "Rcpp::compileAttributes()"
```

## Architecture

### Core Algorithm Flow

1.  **Canonicalize input** (Z/t/log10p → standardized Z-equivalent
    field)
2.  **Build parcel roots** from atlas labels with prior mass caching
3.  **Score parcels** using prior-weighted statistics:
    - κ=0: Variance-stabilized diffuse score
      `U_0(R) = Σπ(v)Z(v) / √Σπ(v)²`
    - κ\>0: Soft-max score `S_κ(R) = (1/κ)·log(Σπ(v)exp(κZ(v)) / Σπ(v))`
    - Combined: `T(R) = max(U_0, max_κ S_κ)`
4.  **Westfall-Young step-down** at parcel level for FWER control
5.  **Recursive octree descent** into significant parcels with
    α-spending
6.  **One-pass sibling scoring** (Rcpp) for children at each split

### Key Statistics

| Statistic | Formula                               | Use Case                      |
|-----------|---------------------------------------|-------------------------------|
| `U_0(R)`  | `Σπ(v)Z(v) / √Σπ(v)²`                 | Distributed/diffuse signals   |
| `S_κ(R)`  | `(1/κ)·log(Σπ(v)·exp(κZ(v)) / Σπ(v))` | Focal signals (κ large → max) |
| `T(R)`    | `max(U_0, max_κ S_κ)`                 | Omnibus score for any pattern |

### Directory Structure

    neurothresh/
    ├── R/                      # R source files
    │   ├── canonicalize.R      # Z/t/log10p → Zeq conversion
    │   ├── hier_scan.R         # Main entry point
    │   ├── score_node.R        # Node scoring functions
    │   ├── wy_stepdown.R       # Westfall-Young step-down
    │   ├── rft_*.R             # RFT baseline methods
    │   ├── tfce.R              # TFCE implementation
    │   └── cluster_fdr.R       # Cluster-FDR methods
    ├── src/                    # Rcpp C++ code
    │   ├── octree_split_info.cpp
    │   ├── score_children_onepass.cpp
    │   ├── split_indices_bbox.cpp
    │   └── RcppExports.cpp
    ├── tests/testthat/         # testthat tests
    └── docs/discussion/        # Design documents (parts 1-5)

### Rcpp C++ Functions

Three critical Rcpp exports for performance: -
`octree_split_info_cpp()` - Build child assignment and cache π/π²
constants - `score_children_onepass_cpp()` - Score all ≤8 children in
single pass (streaming log-sum-exp) - `split_indices_bbox_mass_cpp()` -
Split indices + compute bboxes + π masses in one pass

### neuroim2 Integration Points

``` r
# Key classes from neuroim2
DenseNeuroVol      # 3D statistic maps
LogicalNeuroVol    # Brain masks
ClusteredNeuroVol  # Parcellations (from neuroatlas)
NeuroSpace         # Coordinate system

# Key functions
index_to_grid(vol, idx)  # Linear index → (i,j,k)
grid_to_index(vol, coords)  # (i,j,k) → linear index
conn_comp(vol, threshold)  # Connected components for clustering
```

### neuroatlas Integration

``` r
# Load parcellation (returns ClusteredNeuroVol)
atlas <- neuroatlas::get_schaefer_atlas(parcels="200", networks="7")
cvol <- atlas$atlas

# Access parcel membership
lab_vec <- as.integer(atlas_labels[mask_idx])  # Parcel IDs in mask-space
```

## Algorithm Design Notes

### Why LR-MFT Beats RFT Peak-Height

RFT peak-height is a hard-max detector optimal for “one dominant peak.”
LR-MFT’s soft evidence aggregation is NP-optimal for multi-focal
alternatives (multiple activation blobs), which is biologically
realistic for task fMRI.

### Variance Stabilization (U_0)

Without scaling, singleton voxels have unit variance while large regions
have tiny variance under the null. The `√n_eff` scaling
(`n_eff = Π(R)²/Π_2(R)`) ensures fair comparison across scales.

### Hierarchical α-Spending

Budget allocation: `α_child = (α_node - α_test) · w_child` where
`w_child ∝ π(child)`. This gives a clean FWER proof via union bound +
nesting.

### Computational Optimization

The one-pass sibling scoring eliminates O(m) vector allocations per
permutation by streaming through parent voxels and updating ≤8 child
accumulators simultaneously. This is the key to making permutation-based
inference tractable.

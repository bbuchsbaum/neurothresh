# PRD: neurothresh - Neuroimaging Statistical Thresholding Package

## Introduction

`neurothresh` is an R package for statistical thresholding of neuroimaging data, implementing **LR-MFT (Likelihood-Ratio Matched-Filter Thresholding)** with hierarchical stepdown inference. The package provides a theoretically superior alternative to traditional GRF/RFT methods, along with baseline implementations (RFT, TFCE, cluster-FDR) for comparison.

**Primary Use Case:** Group-level fMRI research analysis (SPM-style workflows)

**Target Users:** Neuroimaging experts comfortable with R and statistics

**Key Innovation:** LR-MFT is provably more powerful than GRF for detecting spatially extended activation blobs, using Neyman-Pearson optimal likelihood-ratio tests with prior-weighted evidence aggregation.

---

## Goals

- Implement LR-MFT with prior-weighted soft evidence statistic T_κ(R)
- Provide hierarchical (top-down) inference with α-spending for multi-scale detection
- Support Westfall-Young step-down for strong FWER control with high power
- Include baseline methods (RFT, TFCE, cluster-FDR) for benchmarking
- Integrate optionally with neuroim2/neuroatlas for parcellation-based priors
- Support flat (uniform) priors for voxel-only analyses
- Provide parametric approximations for speed (permutation optional)
- Handle multiple input stat types: Z-scores, t-statistics, -log10(p)

---

## User Stories

### US-001: Package Infrastructure Setup
**Description:** As a developer, I need the R package skeleton with Rcpp support so that C++ code compiles correctly.

**Acceptance Criteria:**
- [ ] DESCRIPTION file with proper metadata, Imports (Rcpp, neuroim2), Suggests (neuroatlas)
- [ ] NAMESPACE with useDynLib and Rcpp exports
- [ ] src/Makevars configured for C++11
- [ ] src/RcppExports.cpp generated
- [ ] R/RcppExports.R generated
- [ ] `devtools::check()` passes with no errors

---

### US-002: Stat Map Canonicalization
**Description:** As a user, I want to input Z-scores, t-statistics, or -log10(p) values so that the package handles my data format.

**Acceptance Criteria:**
- [ ] `canonicalize_stat(vol, type, df=NULL)` function implemented
- [ ] Supports type = "Z", "t", "neglog10p"
- [ ] t-statistics converted to Z via `qnorm(pt(t, df))`
- [ ] -log10(p) converted to Z via `qnorm(1 - 10^(-x))`
- [ ] Returns NeuroVol with Z-scores
- [ ] Roxygen2 documentation complete
- [ ] Unit tests cover all three input types

---

### US-003: Core C++ Bbox Helpers
**Description:** As a developer, I need C++ functions for bounding box operations so that octree splitting is efficient.

**Acceptance Criteria:**
- [ ] `bbox_from_indices_cpp(indices, dims)` returns 6-element IntegerVector (i0,i1,j0,j1,k0,k1)
- [ ] `split_indices_bbox_cpp(indices, dims)` returns List of 8 child index vectors
- [ ] `split_indices_bbox_mass_cpp(indices, dims, prior_vec)` returns List with indices + prior masses
- [ ] All functions registered in RcppExports
- [ ] Unit tests verify correct octant assignment

---

### US-004: One-Pass Sibling Scoring (Rcpp)
**Description:** As a developer, I need the optimized one-pass scoring algorithm so that all 8 children are scored in a single pass through parent voxels.

**Acceptance Criteria:**
- [ ] `octree_split_info_cpp(indices, dims)` returns octant assignments + bbox
- [ ] `score_children_onepass_cpp(indices, octants, z_vals, prior_vals, kappa)` returns 8 scores
- [ ] Single O(n) pass through parent indices
- [ ] Scores computed as log(sum(pi * exp(kappa * z)))
- [ ] Unit tests verify scores match naive per-child computation

---

### US-005: Prior-Weighted Statistic T_κ(R)
**Description:** As a user, I want to compute the prior-weighted soft evidence statistic for any region R.

**Acceptance Criteria:**
- [ ] `score_set(indices, z_vol, prior_vol, kappa)` R function implemented
- [ ] Uses C++ backend for performance
- [ ] Formula: T_κ(R) = log Σ_{v∈R} π(v) · exp(κ · Z(v))
- [ ] Handles edge cases (empty set, all-zero priors)
- [ ] Unit tests with known values

---

### US-006: Variance-Stabilized Score U_0(R)
**Description:** As a user, I want a fair comparison across regions of different sizes using the variance-stabilized statistic.

**Acceptance Criteria:**
- [ ] `score_set_stabilized(indices, z_vol, prior_vol)` implemented
- [ ] Formula: U_0(R) = Σ π(v)·Z(v) / sqrt(Σ π(v)²)
- [ ] n_eff = (Σπ)² / Σπ² computed for reference
- [ ] Unit tests verify stabilization works across set sizes

---

### US-007: Octree Split Logic
**Description:** As a developer, I need octree splitting to recursively partition voxel sets into dyadic cubes.

**Acceptance Criteria:**
- [ ] `octree_split(indices, dims, min_voxels=8)` R function
- [ ] Splits bbox into 8 octants at midpoints
- [ ] Returns list of child index vectors (non-empty only)
- [ ] Stops splitting when region < min_voxels
- [ ] Uses C++ `split_indices_bbox_cpp` backend

---

### US-008: Permutation Null Distribution
**Description:** As a user, I want to generate null distributions via sign-flipping for FWER control.

**Acceptance Criteria:**
- [ ] `generate_null_scores(z_vol, prior_vol, regions, n_perm, kappa)` implemented
- [ ] Sign-flipping: multiply Z by random ±1 per permutation
- [ ] Returns matrix: n_perm × n_regions
- [ ] Optional: parametric approximation using Gaussian tail bounds
- [ ] Seed parameter for reproducibility

---

### US-009: Westfall-Young Step-Down
**Description:** As a user, I want high-power FWER control via the step-down maxT procedure.

**Acceptance Criteria:**
- [ ] `wy_stepdown(observed_scores, null_matrix, alpha=0.05)` implemented
- [ ] Sort regions by observed score (descending)
- [ ] Compute successive maxima from null matrix
- [ ] Adjust p-values to be monotonic
- [ ] Returns data.frame with region, score, p_adj, rejected
- [ ] Unit tests verify FWER control

---

### US-010: Hierarchical α-Spending
**Description:** As a user, I want top-down hierarchical inference that allocates α budget across tree levels.

**Acceptance Criteria:**
- [ ] `hier_test(node, alpha, null_dist)` recursive function
- [ ] Only tests children if parent rejected
- [ ] Bonferroni split: α_child = α_parent / n_children
- [ ] Optional: Holm or Hochberg upgrade
- [ ] Returns tree structure with rejection status per node
- [ ] Unit tests verify α-spending is valid

---

### US-011: Main Entry Point - hier_scan()
**Description:** As a user, I want a single function to run the full hierarchical LR-MFT analysis.

**Acceptance Criteria:**
- [ ] `hier_scan(z_vol, prior_vol=NULL, parcels=NULL, alpha=0.05, kappa=1, n_perm=1000, method="stepdown")` implemented
- [ ] If prior_vol=NULL, uses flat priors (uniform)
- [ ] If parcels provided, initializes tree from parcel roots
- [ ] Otherwise, starts from whole-brain and uses octree
- [ ] Returns S3 object of class "neurothresh_result" with:
  - significant_regions (list of index vectors)
  - tree structure
  - scores and p-values
  - parameters used
- [ ] Print and summary methods implemented
- [ ] Roxygen2 documentation with examples

---

### US-012: RFT Peak-Level FWER (Baseline)
**Description:** As a user, I want traditional RFT peak-level inference for comparison.

**Acceptance Criteria:**
- [ ] `rft_peak_fwer(z_vol, mask, fwhm, alpha=0.05)` implemented
- [ ] Computes expected Euler characteristic
- [ ] Uses Gaussian field approximation
- [ ] Returns threshold and significant peak coordinates
- [ ] Matches SPM results for validation

---

### US-013: RFT Cluster-Level FWER (Baseline)
**Description:** As a user, I want traditional RFT cluster-level inference for comparison.

**Acceptance Criteria:**
- [ ] `rft_cluster_fwer(z_vol, mask, fwhm, cluster_thresh, alpha=0.05)` implemented
- [ ] Primary threshold at cluster_thresh (e.g., Z>3.1)
- [ ] Cluster size distribution from RFT
- [ ] Returns significant clusters with extent and p-values

---

### US-014: TFCE Transform (Baseline)
**Description:** As a user, I want FSL-style TFCE for threshold-free cluster enhancement.

**Acceptance Criteria:**
- [ ] `tfce_transform(z_vol, mask, H=2, E=0.5, dh=0.1)` implemented
- [ ] Integrates cluster extent over threshold range
- [ ] Formula: TFCE(v) = ∫ e(h)^E · h^H dh
- [ ] Returns transformed NeuroVol
- [ ] Unit tests verify transform properties

---

### US-015: TFCE FWER (Baseline)
**Description:** As a user, I want TFCE with permutation-based FWER control.

**Acceptance Criteria:**
- [ ] `tfce_fwer(z_vol, mask, n_perm=5000, alpha=0.05)` implemented
- [ ] Applies TFCE to observed and permuted data
- [ ] MaxT procedure across permutations
- [ ] Returns threshold and significant voxel mask

---

### US-016: Cluster-FDR (Baseline)
**Description:** As a user, I want topological FDR control across clusters.

**Acceptance Criteria:**
- [ ] `cluster_fdr(z_vol, mask, cluster_thresh, q=0.05, p_method="rft")` implemented
- [ ] Forms clusters at primary threshold
- [ ] Computes cluster p-values (RFT or permutation)
- [ ] Applies BH procedure to cluster p-values
- [ ] Returns significant clusters

---

### US-017: Two-Sided Inference Wrappers
**Description:** As a user, I want to detect both activations and deactivations.

**Acceptance Criteria:**
- [ ] `hier_scan_twosided(z_vol, ...)` runs positive and negative tails
- [ ] Bonferroni correction: α/2 per tail
- [ ] Returns combined results with direction labels
- [ ] All baseline methods also have _twosided variants

---

### US-018: neuroim2 Integration
**Description:** As a user, I want seamless integration with neuroim2 volume classes.
**Acceptance Criteria:**
- [ ] All functions accept NeuroVol objects
- [ ] Mask extraction via LogicalNeuroVol
- [ ] Index ↔ coordinate conversion via NeuroSpace
- [ ] Results returned as NeuroVol where appropriate
- [ ] ClusteredNeuroVol supported for parcellation input

---

### US-019: neuroatlas Integration (Optional)
**Description:** As a user, I want to optionally use atlas parcellations as region priors.

**Acceptance Criteria:**
- [ ] `load_atlas_priors(atlas)` extracts ClusteredNeuroVol
- [ ] Converts parcels to prior weights (e.g., uniform within parcel)
- [ ] `hier_scan(..., parcels=atlas$atlas)` uses parcels as tree roots
- [ ] Works without neuroatlas installed (graceful degradation)

---

### US-020: Result Visualization
**Description:** As a user, I want to visualize thresholding results.

**Acceptance Criteria:**
- [ ] `plot.neurothresh_result(x, vol, ...)` method
- [ ] Creates masked volume with significant regions
- [ ] Compatible with neuroim2 plotting functions
- [ ] Optional: overlay on template

---

### US-021: Comprehensive Test Suite
**Description:** As a developer, I need comprehensive tests for all functionality.

**Acceptance Criteria:**
- [ ] testthat tests for all exported functions
- [ ] Test fixtures with synthetic data
- [ ] Edge case coverage (empty regions, single voxel, etc.)
- [ ] Numerical accuracy tests against known values
- [ ] `devtools::test()` passes with >80% coverage

---

### US-022: Documentation and Examples
**Description:** As a user, I need clear documentation to use the package.

**Acceptance Criteria:**
- [ ] All exported functions have roxygen2 documentation
- [ ] Package-level documentation (?neurothresh)
- [ ] At least one vignette showing typical workflow
- [ ] Examples in function documentation are runnable

---

## Functional Requirements

- **FR-01:** Package compiles with Rcpp on Linux, macOS, Windows
- **FR-02:** All C++ code uses Rcpp best practices (no raw pointers, proper SEXP handling)
- **FR-03:** `canonicalize_stat()` converts Z/t/neglog10p to Z-scores
- **FR-04:** `score_set()` computes T_κ(R) = log Σ π·exp(κ·Z)
- **FR-05:** `score_set_stabilized()` computes U_0(R) with n_eff normalization
- **FR-06:** `octree_split()` partitions indices into 8 octants recursively
- **FR-07:** One-pass scoring scores all 8 children in O(n) time
- **FR-08:** `generate_null_scores()` produces null distribution via sign-flipping
- **FR-09:** `wy_stepdown()` implements Westfall-Young with monotonic p-values
- **FR-10:** `hier_test()` implements α-spending with parent-gated child tests
- **FR-11:** `hier_scan()` is the main entry point combining all components
- **FR-12:** `rft_peak_fwer()` and `rft_cluster_fwer()` implement GRF inference
- **FR-13:** `tfce_transform()` and `tfce_fwer()` implement FSL-style TFCE
- **FR-14:** `cluster_fdr()` implements topological FDR with BH procedure
- **FR-15:** Two-sided wrappers apply Bonferroni α/2 correction per tail
- **FR-16:** All functions accept neuroim2 NeuroVol objects
- **FR-17:** Parcellation priors from neuroatlas are optional, not required
- **FR-18:** Flat (uniform) priors used when no prior_vol provided

---

## Non-Goals (Out of Scope for v1.0)

- Real-time or single-subject inference optimization
- GUI or Shiny interface
- Integration with Python (nibabel, nilearn)
- Automatic FWHM estimation from data
- Surface-based (cortical) analysis
- Bayesian model selection or posterior inference
- Multi-subject random effects modeling (assume input is group stat map)
- HPC/cluster computing parallelization (single-machine only)

---

## Technical Considerations

### Dependencies
- **Required:** Rcpp (≥1.0.0), neuroim2 (≥0.8.0)
- **Suggested:** neuroatlas (≥0.1.0), testthat, knitr, rmarkdown

### Performance
- C++ implementations for all inner loops
- One-pass sibling scoring avoids 8x redundant iterations
- Bounding box operations in C++ avoid R overhead
- Target: 10,000 permutations on 100k voxels in <60 seconds

### Code Organization
```
R/
  neurothresh-package.R    # Package documentation
  canonicalize.R           # Stat type conversion
  scoring.R                # T_κ, U_0 computations
  octree.R                 # Octree splitting logic
  permutation.R            # Null distribution generation
  stepdown.R               # Westfall-Young procedure
  hierarchical.R           # α-spending tree traversal
  hier_scan.R              # Main entry point
  rft.R                    # RFT baseline methods
  tfce.R                   # TFCE baseline methods
  cluster_fdr.R            # Cluster-FDR baseline
  twosided.R               # Two-sided wrappers
  utils.R                  # Helper functions

src/
  bbox_helpers.cpp         # Bounding box operations
  split_indices.cpp        # Octree splitting
  score_children.cpp       # One-pass sibling scoring
  tfce.cpp                 # TFCE transform (optional)
  RcppExports.cpp          # Auto-generated
```

### Reference Materials
Implementation details are in:
- `docs/discussion/part-01-grf-vs-lrmft.md` - Theory and proofs
- `docs/discussion/part-02-set-indexed-inference.md` - T_κ definition, pseudocode
- `docs/discussion/part-03-hierarchical-stepdown.md` - R pseudocode, Westfall-Young
- `docs/discussion/part-04-onepass-rcpp.md` - **Production C++ code (USE VERBATIM)**
- `docs/discussion/part-05-bbox-helpers-and-baselines.md` - **Production C++ code (USE VERBATIM)**
- `lib_notes.md` - neuroim2/neuroatlas API reference

---

## Success Metrics

- `R CMD check --as-cran` passes with 0 errors, 0 warnings
- All testthat tests pass (>80% coverage)
- `hier_scan()` produces valid FWER-controlled results
- LR-MFT shows improved power over RFT in simulation (per discussion notes)
- Package installable via `devtools::install_github()`

---

## Open Questions

1. Should we include simulation functions for power analysis?
2. What default κ value provides best power/robustness tradeoff?
3. Should permutation be truly optional or strongly recommended?
4. Include connectivity constraint for octree (26-connected)?

---

## Implementation Priority

**Phase 1 - Core Infrastructure:**
US-001, US-002, US-003, US-004

**Phase 2 - Core Algorithm:**
US-005, US-006, US-007, US-008, US-009, US-010, US-011

**Phase 3 - Baselines:**
US-012, US-013, US-014, US-015, US-016

**Phase 4 - Integration & Polish:**
US-017, US-018, US-019, US-020, US-021, US-022

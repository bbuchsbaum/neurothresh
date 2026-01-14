# neurothresh Development Task List

**IMPORTANT:** Mark each task [x] IMMEDIATELY after completing it. Do not batch updates.

---

## Phase 1: Infrastructure (DONE - verify and mark)

- [x] **P1-01:** Verify DESCRIPTION, NAMESPACE, src/Makevars exist and are correct
- [x] **P1-02:** Run `Rscript -e "Rcpp::compileAttributes()"` and verify src/RcppExports.cpp generated
- [x] **P1-03:** Run `Rscript -e "devtools::load_all()"` and verify package loads

## Phase 2: Core Functions (verify existing, add docs)

- [x] **P2-01:** Verify `canonicalize_stat()` exists in R/canonicalize.R, add roxygen2 docs, run `devtools::document()`
- [x] **P2-02:** Verify `score_set()` exists in R/scoring.R, add roxygen2 docs
- [x] **P2-03:** Verify `score_set_stabilized()` exists in R/scoring.R, add roxygen2 docs
- [x] **P2-04:** Verify `octree_split()` exists in R/octree.R, add roxygen2 docs
- [x] **P2-05:** Verify `generate_null_scores()` exists in R/permutation.R, add roxygen2 docs
- [x] **P2-06:** Verify `wy_stepdown()` exists in R/stepdown.R, add roxygen2 docs
- [x] **P2-07:** Verify `hier_scan()` exists in R/hier_scan.R, has roxygen2 docs

## Phase 3: Baseline Methods (verify existing, add docs)

- [x] **P3-01:** Verify `rft_peak_fwer()` exists in R/rft.R, add roxygen2 docs
- [x] **P3-02:** Verify `rft_cluster_fwer()` exists in R/rft.R, add roxygen2 docs
- [x] **P3-03:** Verify `tfce_transform()` exists in R/tfce.R, add roxygen2 docs
- [x] **P3-04:** Verify `tfce_fwer()` exists in R/tfce.R, add roxygen2 docs
- [x] **P3-05:** Verify `cluster_fdr()` exists in R/cluster_fdr.R, add roxygen2 docs

## Phase 4: Documentation Generation

- [x] **P4-01:** Run `devtools::document()` to generate all man/*.Rd files
- [x] **P4-02:** Run `R CMD check . --no-manual` and fix any warnings about undocumented objects
- [x] **P4-03:** Ensure all exported functions have @export tag

## Phase 5: Final Verification

- [x] **P5-01:** Run `devtools::test()` - all tests must pass
- [x] **P5-02:** Run `R CMD check . --no-manual --no-vignettes` - must have 0 errors, 0 warnings
- [x] **P5-03:** Run `./verify_completion.sh` - must output "VERIFICATION PASSED"

## COMPLETION GATE

- [x] **GATE:** All above tasks [x] AND verify_completion.sh passes

---

## Notes

- Most code already exists from previous implementation
- Focus is on DOCUMENTATION - adding roxygen2 comments and generating .Rd files
- R CMD check warnings about "Undocumented code objects" must be fixed
- After each task, IMMEDIATELY mark it [x] before moving on

## Ad Hoc Reviews

- [x] **R-001:** Review codebase completeness vs docs/discussion and library notes
- [x] **R-002:** Implement prior normalization and hierarchical alpha-spending mode
- [x] **R-003:** Fix kappa handling, voxel spacing usage, TFCE outputs, and add tests
- [x] **R-004:** Add vignette for neurothresh overview and usage
- [x] **R-005:** Add simulation vignette and sanity-check tests for FWER/power

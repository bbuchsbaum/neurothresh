# neurothresh Development Task List

**IMPORTANT:** Mark each task \[x\] IMMEDIATELY after completing it. Do
not batch updates.

------------------------------------------------------------------------

## Phase 1: Infrastructure (DONE - verify and mark)

**P1-01:** Verify DESCRIPTION, NAMESPACE, src/Makevars exist and are
correct

**P1-02:** Run `Rscript -e "Rcpp::compileAttributes()"` and verify
src/RcppExports.cpp generated

**P1-03:** Run `Rscript -e "devtools::load_all()"` and verify package
loads

## Phase 2: Core Functions (verify existing, add docs)

**P2-01:** Verify
[`canonicalize_stat()`](https://bbuchsbaum.github.io/neurothresh/reference/canonicalize_stat.md)
exists in R/canonicalize.R, add roxygen2 docs, run
`devtools::document()`

**P2-02:** Verify
[`score_set()`](https://bbuchsbaum.github.io/neurothresh/reference/score_set.md)
exists in R/scoring.R, add roxygen2 docs

**P2-03:** Verify
[`score_set_stabilized()`](https://bbuchsbaum.github.io/neurothresh/reference/score_set_stabilized.md)
exists in R/scoring.R, add roxygen2 docs

**P2-04:** Verify
[`octree_split()`](https://bbuchsbaum.github.io/neurothresh/reference/octree_split.md)
exists in R/octree.R, add roxygen2 docs

**P2-05:** Verify
[`generate_null_scores()`](https://bbuchsbaum.github.io/neurothresh/reference/generate_null_scores.md)
exists in R/permutation.R, add roxygen2 docs

**P2-06:** Verify
[`wy_stepdown()`](https://bbuchsbaum.github.io/neurothresh/reference/wy_stepdown.md)
exists in R/stepdown.R, add roxygen2 docs

**P2-07:** Verify
[`hier_scan()`](https://bbuchsbaum.github.io/neurothresh/reference/hier_scan.md)
exists in R/hier_scan.R, has roxygen2 docs

## Phase 3: Baseline Methods (verify existing, add docs)

**P3-01:** Verify
[`rft_peak_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_peak_fwer.md)
exists in R/rft.R, add roxygen2 docs

**P3-02:** Verify
[`rft_cluster_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/rft_cluster_fwer.md)
exists in R/rft.R, add roxygen2 docs

**P3-03:** Verify
[`tfce_transform()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_transform.md)
exists in R/tfce.R, add roxygen2 docs

**P3-04:** Verify
[`tfce_fwer()`](https://bbuchsbaum.github.io/neurothresh/reference/tfce_fwer.md)
exists in R/tfce.R, add roxygen2 docs

**P3-05:** Verify
[`cluster_fdr()`](https://bbuchsbaum.github.io/neurothresh/reference/cluster_fdr.md)
exists in R/cluster_fdr.R, add roxygen2 docs

## Phase 4: Documentation Generation

**P4-01:** Run `devtools::document()` to generate all man/\*.Rd files

**P4-02:** Run `R CMD check . --no-manual` and fix any warnings about
undocumented objects

**P4-03:** Ensure all exported functions have @export tag

## Phase 5: Final Verification

**P5-01:** Run `devtools::test()` - all tests must pass

**P5-02:** Run `R CMD check . --no-manual --no-vignettes` - must have 0
errors, 0 warnings

**P5-03:** Run `./verify_completion.sh` - must output “VERIFICATION
PASSED”

## COMPLETION GATE

**GATE:** All above tasks \[x\] AND verify_completion.sh passes

------------------------------------------------------------------------

## Notes

- Most code already exists from previous implementation
- Focus is on DOCUMENTATION - adding roxygen2 comments and generating
  .Rd files
- R CMD check warnings about “Undocumented code objects” must be fixed
- After each task, IMMEDIATELY mark it \[x\] before moving on

## Ad Hoc Reviews

**R-001:** Review codebase completeness vs docs/discussion and library
notes

**R-002:** Implement prior normalization and hierarchical alpha-spending
mode

**R-003:** Fix kappa handling, voxel spacing usage, TFCE outputs, and
add tests

**R-004:** Add vignette for neurothresh overview and usage

**R-005:** Add simulation vignette and sanity-check tests for FWER/power

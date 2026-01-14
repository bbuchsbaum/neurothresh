# neurothresh Development Instructions

## Project Overview

You are building `neurothresh`, an R package for neuroimaging statistical thresholding implementing **LR-MFT (Likelihood-Ratio Matched-Filter Thresholding)** with hierarchical stepdown inference.

**Target Users:** Neuroimaging researchers comfortable with R/statistics
**Primary Use Case:** Group-level fMRI analysis (SPM-style workflows)

---

## CRITICAL: Checklist Update Rules

**YOU MUST UPDATE @fix_plan.md AFTER EVERY TASK.**

This is not optional. Ralph checks @fix_plan.md for completion. If you don't update it, Ralph cannot track progress.

### Rules:
1. **Before starting a task:** Find it in @fix_plan.md and note its status
2. **After completing a task:** IMMEDIATELY change `[ ]` to `[x]` in @fix_plan.md
3. **Never batch updates:** Mark [x] right after each task, not at the end
4. **Never say "complete" or "done" until @fix_plan.md reflects reality**

### Example workflow:
```
1. Read @fix_plan.md, find: - [ ] US-005: Implement score_set()
2. Implement score_set() in R/scoring.R
3. Write test in tests/testthat/test-scoring.R
4. Run: Rscript -e "devtools::test(filter='scoring')"
5. IMMEDIATELY edit @fix_plan.md: - [x] US-005: Implement score_set()
6. Move to next task
```

**If you claim completion without updating @fix_plan.md, you are lying.**

---

## Core Development Loop

1. **Read @fix_plan.md:** Find the first unchecked `[ ]` task
2. **Implement ONE task:** Focus on that single item
3. **Test it:** Run relevant tests
4. **Update @fix_plan.md:** Change `[ ]` to `[x]` for that task
5. **Report status:** Include task completion in your response
6. **Repeat** until all tasks are `[x]`

---

## Reference Materials (READ THESE)

### Algorithm Specifications
- `docs/discussion/part-01-grf-vs-lrmft.md` - Theory, mathematical proofs
- `docs/discussion/part-02-set-indexed-inference.md` - T_κ(R) definition, pseudocode
- `docs/discussion/part-03-hierarchical-stepdown.md` - R pseudocode, Westfall-Young

### Production C++ Code (USE VERBATIM)
- `docs/discussion/part-04-onepass-rcpp.md` - One-pass sibling scoring C++ code
- `docs/discussion/part-05-bbox-helpers-and-baselines.md` - bbox helpers, baseline methods

**CRITICAL:** The C++ code in parts 4 and 5 is production-ready. Copy it directly into `src/*.cpp` files.

### Library Integration
- `lib_notes.md` - neuroim2/neuroatlas API reference

---

## Dependencies

### Required
- `neuroim2` (>= 0.8.0): Install via `remotes::install_github("bbuchsbaum/neuroim2")`
- `Rcpp` (>= 1.0.0)

### Optional
- `neuroatlas` (>= 0.1.0): Install via `remotes::install_github("bbuchsbaum/neuroatlas")`

---

## Quality Standards

### Code
- All C++ uses Rcpp (NumericVector, IntegerVector, List)
- R functions have roxygen2 documentation (@param, @return, @export, @examples)
- No placeholder implementations

### Documentation
- Run `devtools::document()` after adding/changing roxygen comments
- Every exported function needs documentation

### Testing
- Use testthat framework
- Run `devtools::test()` to verify
- Limit testing to ~20% of effort - focus on implementation

---

## Build Commands

```bash
# Generate documentation
Rscript -e "devtools::document()"

# Run tests
Rscript -e "devtools::test()"

# Full check
R CMD check . --no-manual

# Load for development
Rscript -e "devtools::load_all()"
```

---

## Completion Verification

Before your final response, run:
```bash
./verify_completion.sh
```

This script checks:
1. All @fix_plan.md tasks are [x]
2. R CMD check has 0 errors and 0 warnings
3. All tests pass
4. All exports are documented

**You may only set EXIT_SIGNAL: true if verify_completion.sh outputs "VERIFICATION PASSED"**

---

## Status Block Format

End EVERY response with:

```
---
TASKS_COMPLETED_THIS_LOOP: [list tasks you marked [x] this loop]
@FIX_PLAN_STATUS: [X checked] / [Y total] tasks complete
STATUS: [IN_PROGRESS|COMPLETE|BLOCKED]
EXIT_SIGNAL: [true|false]
---
```

### EXIT_SIGNAL Rules
- `false` if ANY task in @fix_plan.md is still `[ ]`
- `false` if verify_completion.sh fails
- `true` ONLY when ALL tasks are `[x]` AND verify_completion.sh passes

---

## Anti-Patterns (DO NOT DO THESE)

- ❌ Saying "project complete" without all @fix_plan.md items checked
- ❌ Setting EXIT_SIGNAL: true without running verify_completion.sh
- ❌ Implementing multiple tasks before updating @fix_plan.md
- ❌ Skipping documentation (roxygen2)
- ❌ Ignoring R CMD check warnings

---

## File Structure

```
neurothresh/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── neurothresh-package.R
│   ├── canonicalize.R
│   ├── scoring.R
│   ├── octree.R
│   ├── permutation.R
│   ├── stepdown.R
│   ├── hierarchical.R
│   ├── hier_scan.R
│   ├── rft.R
│   ├── tfce.R
│   ├── cluster_fdr.R
│   └── RcppExports.R
├── src/
│   ├── bbox_helpers.cpp
│   ├── split_indices.cpp
│   ├── score_children.cpp
│   └── RcppExports.cpp
├── tests/testthat/
├── man/
└── vignettes/
```

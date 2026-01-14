# Neurothresh Implementation Plan with Ralph

## Overview

We will use the **Ralph autonomous development framework** to implement the `neurothresh` R package end-to-end. Ralph operates in iterative loops, reading from structured files and autonomously implementing until completion.

---

## Phase 1: Create the PRD

**Goal**: Transform our captured discussion (parts 1-5) into a structured Product Requirements Document.

**Action**: Run `/prd` skill to generate a comprehensive PRD that captures:
- Core algorithms (LR-MFT, hierarchical stepdown, Westfall-Young)
- C++/Rcpp implementation requirements
- Integration with neuroim2 and neuroatlas
- Standard baseline methods (RFT, TFCE, cluster-FDR)
- API surface and user-facing functions

**Output**: `docs/PRD.md`

---

## Phase 2: Initialize Ralph Project Structure

**Goal**: Set up Ralph's required file structure for autonomous development.

### Files to Create

```
neurothresh/
├── PROMPT.md           # Main development instructions for Ralph
├── @fix_plan.md        # Prioritized task list (Ralph-controlled)
├── @AGENT.md           # R package build/test instructions
├── specs/              # Technical specifications
│   ├── algorithms.md   # Core algorithm specs from discussion
│   ├── api.md          # Public API specification
│   └── rcpp.md         # C++ implementation details
├── R/                  # R source files
├── src/                # C++ source files (Rcpp)
├── tests/              # testthat tests
├── man/                # Documentation (roxygen2)
├── DESCRIPTION         # R package metadata
├── NAMESPACE           # Exports
└── docs/
    ├── PRD.md
    └── discussion/     # (already exists)
```

---

## Phase 3: Structure the @fix_plan.md

Ralph reads this file to determine what to work on. Tasks should be organized by priority:

### High Priority (Core Infrastructure)
1. [ ] Set up R package skeleton (DESCRIPTION, NAMESPACE)
2. [ ] Create Rcpp infrastructure (src/, Makevars)
3. [ ] Implement canonicalization layer (Z/t/log10p conversion)
4. [ ] Implement core C++ helpers (bbox, split_indices, score_children)

### High Priority (Core Algorithms)
5. [ ] Implement prior-weighted statistic T_κ(R)
6. [ ] Implement variance-stabilized score U_0(R)
7. [ ] Implement one-pass sibling scoring (Rcpp)
8. [ ] Implement octree split logic

### Medium Priority (Inference Framework)
9. [ ] Implement permutation null distribution
10. [ ] Implement Westfall-Young step-down
11. [ ] Implement hierarchical α-spending
12. [ ] Main entry point: hier_scan()

### Medium Priority (Baseline Methods)
13. [ ] Implement RFT peak/cluster FWER
14. [ ] Implement TFCE transform and FWER
15. [ ] Implement cluster-FDR

### Medium Priority (Integration)
16. [ ] Integration with neuroim2 (NeuroVol, NeuroSpace)
17. [ ] Integration with neuroatlas (parcellations, ROIs)
18. [ ] Two-sided inference wrappers

### Low Priority (Polish)
19. [ ] User-facing convenience wrappers
20. [ ] Vignettes and examples
21. [ ] Performance optimization
22. [ ] Comprehensive test coverage

---

## Phase 4: Configure @AGENT.md for R

Ralph needs R-specific build/test instructions:

```markdown
# Build Instructions

## Setup
```bash
# Install dependencies
Rscript -e "install.packages(c('Rcpp', 'testthat', 'roxygen2', 'devtools'))"
Rscript -e "devtools::install_deps()"
```

## Build
```bash
Rscript -e "devtools::document()"
Rscript -e "devtools::build()"
R CMD INSTALL .
```

## Test
```bash
Rscript -e "devtools::test()"
R CMD check --as-cran .
```

## Development
```bash
Rscript -e "devtools::load_all()"
```
```

---

## Phase 5: Configure PROMPT.md

The main instruction file Ralph reads each loop:

```markdown
# Neurothresh Development Instructions

## Project Context
You are building an R package called `neurothresh` for neuroimaging statistical thresholding.
The package implements LR-MFT (Likelihood-Ratio Matched-Filter Thresholding) with hierarchical
stepdown inference, plus baseline methods (RFT, TFCE, cluster-FDR).

## Dependencies
- neuroim2: ~/code/neuroim2 (NeuroVol, NeuroSpace, indexing)
- neuroatlas: ~/code/neuroatlas (parcellations, ROI definitions)

## Core Loop
1. Study specs/ and PRD.md for requirements
2. Check @fix_plan.md for highest priority incomplete task
3. Implement ONE task per loop
4. Run tests: devtools::test()
5. Update @fix_plan.md marking completed items
6. Provide status block

## Quality Standards
- All C++ code uses Rcpp best practices
- R functions have roxygen2 documentation
- Tests use testthat framework
- No placeholder implementations
- Code must integrate with neuroim2/neuroatlas types

## Status Block Format
STATUS: [IN_PROGRESS|COMPLETE|BLOCKED]
TASK: [current task description]
TESTS: [pass/fail count]
EXIT_SIGNAL: [true|false]
```

---

## Phase 6: Execute Ralph Loop

Once files are in place:

```bash
# Start Ralph autonomous development
ralph --timeout 60 --calls 100 --monitor

# Or without monitoring
ralph --timeout 30
```

Ralph will iterate:
1. Read PROMPT.md
2. Check @fix_plan.md for next task
3. Implement the task
4. Run R CMD check / devtools::test()
5. Update @fix_plan.md
6. Report status
7. Loop until EXIT_SIGNAL=true

---

## Execution Order

| Step | Action | Output |
|------|--------|--------|
| 1 | Run `/prd` skill | `docs/PRD.md` |
| 2 | Create PROMPT.md | Ralph instructions |
| 3 | Create @fix_plan.md | Task checklist |
| 4 | Create @AGENT.md | Build instructions |
| 5 | Create specs/*.md | Algorithm specifications |
| 6 | Create R package skeleton | DESCRIPTION, NAMESPACE |
| 7 | Start Ralph | Autonomous implementation |

---

## Success Criteria

Ralph exits when:
- All @fix_plan.md items marked [x]
- R CMD check passes with no errors/warnings
- All testthat tests pass
- Core functions documented and exported
- Package installable: `devtools::install()`

---

## Notes on neuroim2/neuroatlas Integration

The package will use:
- `NeuroVol` for 3D statistical maps
- `NeuroSpace` for coordinate systems
- `NeuroVec` for vectorized operations
- Parcellation objects from neuroatlas for ROI definitions
- Index conversion utilities for voxel ↔ coordinate mapping

Ralph should read these packages' documentation when implementing integration points.

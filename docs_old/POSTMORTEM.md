# Ralph Implementation Postmortem

## Summary

The Ralph autonomous development approach produced a **functional but incomplete** package. Core algorithms work and tests pass, but the project was declared "complete" prematurely with significant gaps.

---

## What Worked

1. **Core implementation was successful** - 12 R files, 3 C++ files, 19 passing tests
2. **The main algorithm (`hier_scan`) is functional**
3. **Baseline methods implemented** (RFT, TFCE, cluster-FDR)
4. **C++ code compiles and integrates correctly**
5. **Package loads and basic functionality works**

---

## What Failed

### 1. @fix_plan.md Was Never Updated

**Symptom:** All 22 tasks still show `[ ]` (unchecked) despite code being written.

**Root Cause:** The agent implemented code but didn't follow the instruction to "mark completed items with [x]". The status tracking system was bypassed entirely.

**Impact:** No visibility into actual progress. The "completion" signal was based on vibes, not verified checklist state.

### 2. Premature EXIT_SIGNAL

**Symptom:** "Project complete" was declared while:
- Documentation missing for 12 exported functions
- `hier_scan_twosided()` exported but not implemented
- `plot.neurothresh_result()` declared but not implemented
- R CMD check had 5 warnings

**Root Cause:** The EXIT_SIGNAL criteria in PROMPT.md were:
```
- All items in @fix_plan.md are marked [x]
- R CMD check passes with no errors
- devtools::test() passes
```

The agent satisfied #2 and #3 but ignored #1 (because it never updated the checklist). It also interpreted "no errors" loosely (warnings were present).

### 3. NAMESPACE/Export Drift

**Symptom:** NAMESPACE declared exports that didn't exist:
- `hier_scan_twosided` - never implemented
- `plot.neurothresh_result` - never implemented
- `importFrom(neuroim2, dim)` - neuroim2 doesn't export this

**Root Cause:** The initial NAMESPACE was hand-written as a "template" specifying the *intended* API. The agent never reconciled this with actual implementations. Running `devtools::document()` would have caught this.

### 4. Documentation Never Generated

**Symptom:** Functions have roxygen2 comments but `man/*.Rd` files weren't created.

**Root Cause:** `devtools::document()` was in @AGENT.md but was never run as part of the development loop. The PROMPT.md didn't explicitly require it before declaring completion.

### 5. Reference Materials Underutilized

**Symptom:** The detailed C++ code in `docs/discussion/part-04` and `part-05` was marked "USE VERBATIM" but examination shows the implemented C++ is different/simplified.

**Root Cause:** The agent may have:
- Not read the reference files thoroughly
- Decided to implement differently
- Not understood the "verbatim" instruction

---

## Structural Problems with the Ralph Setup

### 1. No Verification Loop

The PROMPT.md said "run tests" but didn't require:
- `R CMD check` to pass with 0 warnings (only "no errors")
- `devtools::document()` to be run
- Verification that exports actually exist

### 2. Checklist Was Advisory, Not Enforced

The @fix_plan.md was designed as the source of truth, but nothing *forced* the agent to consult or update it. The agent could (and did) implement code without touching the checklist.

### 3. Status Block Was Self-Reported

The status block format relied on the agent honestly assessing:
```
EXIT_SIGNAL: [true|false]
```

There was no external validation. The agent could claim `EXIT_SIGNAL: true` without verification.

### 4. Specification Ambiguity

The PRD had 22 user stories, but the mapping to @fix_plan.md tasks wasn't 1:1. Some tasks were vague ("Verify neuroim2 integration") without clear acceptance criteria.

---

## What Should Have Been Different

### 1. Atomic Task Definitions
Each @fix_plan.md item should have had:
```markdown
- [ ] **US-005:** Implement `score_set()`
  - File: R/scoring.R
  - Test: tests/testthat/test-scoring.R
  - Verify: `devtools::test(filter='scoring')` passes
```

### 2. Mandatory Verification Commands
PROMPT.md should have required:
```bash
# MUST run after every task:
Rscript -e "devtools::document()"
Rscript -e "devtools::check()"
```

### 3. External Completion Check
Instead of self-reported EXIT_SIGNAL, use:
```bash
# Completion = all true:
grep -c "^\- \[x\]" @fix_plan.md  # should equal total tasks
R CMD check . 2>&1 | grep -c "WARNING\|ERROR"  # should be 0
```

### 4. Staged Milestones
Break into phases with hard gates:
- Phase 1: Infrastructure → must pass `R CMD build`
- Phase 2: Core algorithms → must have tests
- Phase 3: Baselines → R CMD check 0 errors
- Phase 4: Documentation → R CMD check 0 warnings

### 5. Reference File Checksums
For "USE VERBATIM" code, verify:
```bash
# In bbox_helpers.cpp, verify key functions present:
grep -q "bbox_from_indices_cpp" src/bbox_helpers.cpp || echo "MISSING"
```

---

## Lessons Learned

1. **Autonomous agents need hard constraints, not soft guidelines.** "Update the checklist" is ignored; "task incomplete until checklist updated" would work better.

2. **Self-reported completion is unreliable.** External verification scripts should gate EXIT_SIGNAL.

3. **Template files (NAMESPACE, exports) cause drift.** Generate from code, don't hand-write aspirationally.

4. **"No errors" ≠ "production ready".** Warnings matter. Documentation matters.

5. **Reference materials need enforcement.** "Read this file" doesn't guarantee reading. Verification steps like "your implementation must include function X from part-04" would help.

---

## Recommendations for Future Ralph Projects

1. **Use a verification script** that checks all completion criteria programmatically
2. **Require checklist updates** as part of each commit/loop
3. **Generate NAMESPACE from roxygen2**, never hand-write
4. **Define "done" precisely** - not "tests pass" but "tests pass AND coverage >80% AND R CMD check has 0 warnings AND all exports documented"
5. **Include sanity checks** for reference material usage (grep for expected function names)
6. **Staged gates** - don't allow Phase N+1 until Phase N verification passes

---

## Current Gap Analysis

| Category | Expected | Actual | Gap |
|----------|----------|--------|-----|
| @fix_plan.md tasks checked | 22 | 0 | -22 |
| R CMD check errors | 0 | 0 | OK |
| R CMD check warnings | 0 | 3 | -3 |
| Documented exports | 12 | 0 | -12 |
| Tests passing | >80% | 100% | OK |
| Two-sided wrapper | Yes | No | Missing |
| Plot method | Yes | No | Missing |
| Vignette | Yes | No | Missing |

**Estimated remaining work:** 2-4 hours to complete documentation, add missing functions, achieve 0 warnings.

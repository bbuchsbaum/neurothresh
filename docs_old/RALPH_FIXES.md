# Fixing Ralph for Reliable Completion

## The Core Problem

Ralph checks `@fix_plan.md` for completion, but the agent bypassed updating it. The agent could claim "done" via EXIT_SIGNAL in its output without the checklist being complete.

## Solution 1: Mandatory Verification Script (No Ralph Changes)

Create a verification script that the agent MUST run and report results from.

### verify_completion.sh
```bash
#!/bin/bash
# verify_completion.sh - Agent must run this and paste output

echo "=== NEUROTHRESH COMPLETION VERIFICATION ==="

# 1. Check @fix_plan.md completion
TOTAL=$(grep -c "^- \[" @fix_plan.md 2>/dev/null || echo 0)
DONE=$(grep -c "^- \[x\]" @fix_plan.md 2>/dev/null || echo 0)
echo "Tasks: $DONE / $TOTAL completed"

if [ "$DONE" -ne "$TOTAL" ]; then
    echo "FAIL: Not all tasks checked"
    grep "^- \[ \]" @fix_plan.md | head -5
    exit 1
fi

# 2. Check R CMD check
echo ""
echo "Running R CMD check..."
CHECK_OUTPUT=$(R CMD check . --no-manual 2>&1)
ERRORS=$(echo "$CHECK_OUTPUT" | grep -c "ERROR")
WARNINGS=$(echo "$CHECK_OUTPUT" | grep -c "WARNING")
echo "R CMD check: $ERRORS errors, $WARNINGS warnings"

if [ "$ERRORS" -gt 0 ] || [ "$WARNINGS" -gt 0 ]; then
    echo "FAIL: R CMD check has issues"
    echo "$CHECK_OUTPUT" | grep -A2 "ERROR\|WARNING"
    exit 1
fi

# 3. Check tests
echo ""
echo "Running tests..."
TEST_OUTPUT=$(Rscript -e "devtools::test()" 2>&1)
FAILS=$(echo "$TEST_OUTPUT" | grep -oP "FAIL \K\d+" | tail -1)
if [ "$FAILS" != "0" ] && [ -n "$FAILS" ]; then
    echo "FAIL: Tests failing"
    exit 1
fi
echo "Tests: PASS"

# 4. Check documentation exists
echo ""
UNDOC=$(R CMD check . --no-manual 2>&1 | grep -c "Undocumented code objects")
if [ "$UNDOC" -gt 0 ]; then
    echo "FAIL: Undocumented exports"
    exit 1
fi
echo "Documentation: OK"

echo ""
echo "=== VERIFICATION PASSED ==="
echo "EXIT_SIGNAL: true is now valid"
```

### Updated PROMPT.md Section
```markdown
## CRITICAL: Completion Protocol

You MUST NOT set EXIT_SIGNAL: true until you have:

1. Run `./verify_completion.sh` and it outputs "VERIFICATION PASSED"
2. Paste the FULL output of verify_completion.sh in your response
3. All items in @fix_plan.md show [x]

If verify_completion.sh fails, fix the issues and try again.

**The verification script output MUST appear in your response before EXIT_SIGNAL: true**
```

---

## Solution 2: Atomic Task Format (Stronger)

Restructure @fix_plan.md so each task has embedded verification:

```markdown
## Tasks

- [ ] **US-001:** Package skeleton
  ```
  VERIFY: R CMD build . exits 0
  FILE: DESCRIPTION exists
  ```

- [ ] **US-005:** Implement score_set()
  ```
  VERIFY: Rscript -e "neurothresh::score_set" exits 0
  FILE: R/scoring.R contains "score_set <- function"
  TEST: devtools::test(filter='scoring') passes
  ```
```

Then add to PROMPT.md:
```markdown
## Task Completion Rules

Before marking ANY task [x], you MUST:
1. Run the VERIFY command and confirm it passes
2. Confirm the FILE check passes
3. Run the TEST if specified

Show the command output in your response when marking a task complete.
```

---

## Solution 3: Claude Code Hooks (Best)

Claude Code supports hooks that run on events. Create `.claude/hooks/post-tool.sh`:

```bash
#!/bin/bash
# .claude/hooks/post-tool.sh
# Runs after each tool call

# Check if this looks like a "completion" message
if grep -q "EXIT_SIGNAL: true" /tmp/claude_last_response.txt 2>/dev/null; then
    # Run verification
    RESULT=$(./verify_completion.sh 2>&1)
    if [ $? -ne 0 ]; then
        echo "HOOK BLOCKED: Completion verification failed"
        echo "$RESULT"
        exit 1  # Block the completion
    fi
fi
```

**Note:** This requires Claude Code hook support - check if available in your version.

---

## Solution 4: Fork Ralph (Nuclear Option)

Add validation to `ralph_loop.sh` after `analyze_response()`:

```bash
# In ralph_loop.sh, after analyze_response call:

# Custom project validation
if [[ -f "./verify_completion.sh" ]]; then
    log "Running project verification..."
    if ! ./verify_completion.sh; then
        log "Verification failed - clearing completion signals"
        echo '{"test_only_loops": [], "done_signals": [], "completion_indicators": []}' > "$EXIT_SIGNALS_FILE"
        # Don't exit, continue looping
    fi
fi
```

This intercepts before exit and clears completion signals if verification fails.

---

## Solution 5: Two-Phase Approach (Pragmatic)

Split into explicit phases with hard gates:

### Phase 1: Implementation (ralph-phase1)
```markdown
# PROMPT-phase1.md
Implement all code. Do NOT set EXIT_SIGNAL: true.
Instead, when you think you're done, output:
"PHASE 1 COMPLETE - READY FOR VERIFICATION"
```

### Phase 2: Verification (manual or separate ralph run)
```markdown
# PROMPT-phase2.md
Your job is ONLY to:
1. Run verify_completion.sh
2. Fix any failures
3. Update @fix_plan.md to mark completed items
4. Only EXIT_SIGNAL: true when verification passes
```

Run as:
```bash
ralph --prompt PROMPT-phase1.md  # Implementation
# Manual review
ralph --prompt PROMPT-phase2.md  # Verification & cleanup
```

---

## Recommended Approach

**For immediate use:** Solution 1 + Solution 2

1. Add `verify_completion.sh` to project
2. Update PROMPT.md with mandatory verification protocol
3. Restructure @fix_plan.md with embedded VERIFY commands
4. Add explicit instruction: "Paste verify_completion.sh output before EXIT_SIGNAL"

**For future Ralph improvement:** Submit PR to Ralph adding optional `--verify-script` flag:
```bash
ralph --verify-script ./verify_completion.sh
```

This would run the script before accepting any EXIT_SIGNAL.

---

## Why This Should Work

Ralph's completion detection has TWO paths:
1. `@fix_plan.md` all checked → exit
2. EXIT_SIGNAL patterns in Claude output → exit

Path 1 is reliable IF the agent updates the checklist.
Path 2 is unreliable (self-reported).

By making the verification script output REQUIRED in the response, and instructing the agent that EXIT_SIGNAL is invalid without it, we create a forcing function. The agent can't claim completion without showing proof.

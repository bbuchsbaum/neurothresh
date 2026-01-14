#!/bin/bash
# verify_completion.sh - Must pass before EXIT_SIGNAL: true is valid
# Run this and paste FULL output in your response

echo "╔════════════════════════════════════════════════════════════╗"
echo "║        NEUROTHRESH COMPLETION VERIFICATION                 ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

PASS=0
FAIL=0

# 1. Check @fix_plan.md completion
echo "┌─ 1. Checking @fix_plan.md ─────────────────────────────────┐"
TOTAL=$(grep -c "^- \[" @fix_plan.md 2>/dev/null || echo "0")
DONE=$(grep -c "^- \[x\]" @fix_plan.md 2>/dev/null || echo "0")
echo "   Tasks completed: $DONE / $TOTAL"

if [ "$DONE" = "$TOTAL" ] && [ "$TOTAL" != "0" ]; then
    echo "   ✓ PASS: All tasks checked"
    PASS=$((PASS + 1))
else
    echo "   ✗ FAIL: Unchecked tasks remain:"
    grep "^- \[ \]" @fix_plan.md 2>/dev/null | head -5 | sed 's/^/      /'
    FAIL=$((FAIL + 1))
fi
echo ""

# 2. Check R CMD check
echo "┌─ 2. Running R CMD check ───────────────────────────────────┐"
rm -rf ..Rcheck ../..Rcheck src/*.o src/*.so 2>/dev/null || true
CHECK_OUTPUT=$(LC_ALL=C R CMD check -o .. . --no-manual --no-vignettes 2>&1)
STATUS_LINE=$(echo "$CHECK_OUTPUT" | grep "^Status:")
ERRORS=$(echo "$STATUS_LINE" | grep -o "[0-9]* ERROR" | grep -o "[0-9]*" || echo "0")
WARNINGS=$(echo "$STATUS_LINE" | grep -o "[0-9]* WARNING" | grep -o "[0-9]*" || echo "0")

# Handle empty values
[ -z "$ERRORS" ] && ERRORS="0"
[ -z "$WARNINGS" ] && WARNINGS="0"

echo "   Errors: $ERRORS"
echo "   Warnings: $WARNINGS"

if [ "$ERRORS" = "0" ] && [ "$WARNINGS" = "0" ]; then
    echo "   ✓ PASS: R CMD check clean"
    PASS=$((PASS + 1))
else
    echo "   ✗ FAIL: R CMD check has issues"
    echo "$CHECK_OUTPUT" | grep -A2 "WARNING\|ERROR" | head -15 | sed 's/^/      /'
    FAIL=$((FAIL + 1))
fi
echo ""

# 3. Check tests pass
echo "┌─ 3. Running tests ─────────────────────────────────────────┐"
TEST_OUTPUT=$(Rscript -e "devtools::test()" 2>&1)
# Look for "FAIL X" pattern in results line
TEST_RESULT=$(echo "$TEST_OUTPUT" | grep "FAIL" | tail -1)
if echo "$TEST_RESULT" | grep -q "FAIL 0"; then
    echo "   ✓ PASS: All tests pass"
    PASS=$((PASS + 1))
elif [ -z "$TEST_RESULT" ]; then
    # No FAIL line means all passed
    echo "   ✓ PASS: All tests pass"
    PASS=$((PASS + 1))
else
    echo "   ✗ FAIL: Some tests failing"
    echo "$TEST_RESULT" | sed 's/^/      /'
    FAIL=$((FAIL + 1))
fi
echo ""

# 4. Check exports are documented
echo "┌─ 4. Checking documentation ────────────────────────────────┐"
if echo "$CHECK_OUTPUT" | grep -q "Undocumented code objects"; then
    echo "   ✗ FAIL: Undocumented exports"
    echo "$CHECK_OUTPUT" | grep -A5 "Undocumented code objects" | head -7 | sed 's/^/      /'
    FAIL=$((FAIL + 1))
else
    echo "   ✓ PASS: All exports documented"
    PASS=$((PASS + 1))
fi
echo ""

# 5. Check key functions are loadable
echo "┌─ 5. Checking key exports ──────────────────────────────────┐"
# First install/load the package
Rscript -e "devtools::load_all(quiet=TRUE)" 2>/dev/null

KEY_CHECK=$(Rscript -e "
  devtools::load_all(quiet=TRUE)
  fns <- c('hier_scan', 'score_set', 'canonicalize_stat', 'wy_stepdown')
  missing <- fns[!sapply(fns, exists, where=asNamespace('neurothresh'))]
  if(length(missing) > 0) cat('MISSING:', paste(missing, collapse=' '))
" 2>&1)

if echo "$KEY_CHECK" | grep -q "MISSING:"; then
    echo "   ✗ FAIL: $KEY_CHECK"
    FAIL=$((FAIL + 1))
else
    echo "   ✓ PASS: Key functions exported"
    PASS=$((PASS + 1))
fi
echo ""

# Summary
TOTAL_CHECKS=$((PASS + FAIL))
echo "╔════════════════════════════════════════════════════════════╗"
if [ "$FAIL" -eq 0 ]; then
    echo "║  ✓ VERIFICATION PASSED ($PASS/$TOTAL_CHECKS checks)                       ║"
    echo "║  EXIT_SIGNAL: true is now VALID                            ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    exit 0
else
    echo "║  ✗ VERIFICATION FAILED ($FAIL failures)                        ║"
    echo "║  EXIT_SIGNAL: true is NOT VALID                            ║"
    echo "║  Fix issues above and re-run this script                   ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    exit 1
fi

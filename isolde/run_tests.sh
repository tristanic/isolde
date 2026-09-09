#!/bin/bash
# ===========================================================================
# run_tests.sh -- run every self-contained ISOLDE test under src/tests/
# headlessly (via run_chimerax.sh) and report pass/fail per file.
#
# POSIX counterpart of run_tests.bat, with the same contract: a
# "self-contained" test runs itself when ChimeraX injects `session` (the
# `if session is not None: run(session)` footer) and prints "ALL PASS" on
# success / calls _fail() (which raises SystemExit) otherwise.
#
# Files with no "ALL PASS" sentinel in their source are SKIPPED, not failed --
# e.g. test_simulation.py is the interactive SimTester harness, not a headless
# pass/fail test. Any future test adopting the convention is picked up
# automatically.
#
# Success is judged by "ALL PASS" appearing in the captured output, not by exit
# code: a SystemExit raised inside a ChimeraX --script does not reliably surface
# as a process exit code.
#
# Usage:  ./run_tests.sh [release]
#   release -> forwarded to run_chimerax.sh to target the stable install
#              (omitted -> the daily build), matching run_chimerax.sh's own
#              argument convention.
#
# Optional: ISOLDE_TEST_TIMEOUT=<seconds> caps each test. Unset means no cap,
# matching run_tests.bat. It is implemented by hand because macOS ships no GNU
# `timeout`; note the output must go to a FILE rather than through a pipe, or
# the watchdog's own `sleep` holds the pipe open and the read blocks for the
# full timeout even after the test has exited.
# ===========================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE=""
if [ "${1:-}" = "release" ]; then
    RELEASE="release"
fi

TESTDIR="$HERE/src/tests"
OUTDIR="$(mktemp -d)"
trap 'rm -rf "$OUTDIR"' EXIT

NPASS=0; NFAIL=0; NSKIP=0; FAILED=""

# Run from the bundle dir so the forward-slash --script path resolves the same
# way it does when invoked by hand.
cd "$HERE" || exit 1

echo "==========================================================================="
echo "Running ISOLDE tests in $TESTDIR"
echo "==========================================================================="

run_one() {
    local script="$1" log="$2"
    if [ -n "${ISOLDE_TEST_TIMEOUT:-}" ]; then
        ./run_chimerax.sh $RELEASE --nogui --exit --script "$script" > "$log" 2>&1 &
        local pid=$!
        ( sleep "$ISOLDE_TEST_TIMEOUT"; kill -9 "$pid" 2>/dev/null ) >/dev/null 2>&1 &
        local watcher=$!
        wait "$pid" 2>/dev/null
        # Kill the watchdog AND its sleep child, or the sleep lingers.
        pkill -P "$watcher" 2>/dev/null
        kill -9 "$watcher" 2>/dev/null
        wait "$watcher" 2>/dev/null
    else
        ./run_chimerax.sh $RELEASE --nogui --exit --script "$script" > "$log" 2>&1
    fi
    return 0
}

for f in "$TESTDIR"/test_*.py; do
    [ -e "$f" ] || continue
    name="$(basename "$f")"
    if ! grep -q "ALL PASS" "$f"; then
        echo "[SKIP] $name  (no \"ALL PASS\" sentinel: not a self-running test)"
        NSKIP=$((NSKIP + 1))
        continue
    fi
    echo
    echo "--- $name ---"
    log="$OUTDIR/$name.log"
    run_one "src/tests/$name" "$log"
    cat "$log"
    if grep -q "ALL PASS" "$log"; then
        echo "[PASS] $name"
        NPASS=$((NPASS + 1))
    else
        echo "[FAIL] $name"
        NFAIL=$((NFAIL + 1))
        FAILED="$FAILED $name"
    fi
done

echo
echo "==========================================================================="
echo "Summary: $NPASS passed, $NFAIL failed, $NSKIP skipped"
[ -n "$FAILED" ] && echo "Failed:$FAILED"
echo "==========================================================================="

[ "$NFAIL" -gt 0 ] && exit 1
exit 0

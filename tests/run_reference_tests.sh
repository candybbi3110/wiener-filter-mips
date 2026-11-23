#!/usr/bin/env bash
# Runs the C++ reference model on every test case and checks both
# the output file (expected_output.txt) and the exit status (expected_exit_code).
# Usage: tests/run_reference_tests.sh [path/to/wiener_ref]   (default: build/wiener_ref)
set -u
BIN="${1:-build/wiener_ref}"
DIR="$(cd "$(dirname "$0")" && pwd)"
pass=0; fail=0
for tc in "$DIR"/tc*/; do
  name="$(basename "$tc")"
  out="$(mktemp)"
  "$BIN" "$tc/input.txt" "$tc/desired.txt" "$out" > /dev/null
  code=$?
  want_code="$(cat "$tc/expected_exit_code" 2>/dev/null || echo 0)"
  ok=1
  if [ "$code" -ne "$want_code" ]; then
    echo "FAIL  $name  exit code $code, expected $want_code"; ok=0
  fi
  # whitespace-insensitive compare (MARS output has trailing spaces)
  if ! diff -wB "$out" "$tc/expected_output.txt" > /dev/null; then
    echo "FAIL  $name  output differs:"; diff -wB "$out" "$tc/expected_output.txt"; ok=0
  fi
  if [ "$ok" -eq 1 ]; then echo "PASS  $name"; pass=$((pass+1)); else fail=$((fail+1)); fi
  rm -f "$out"
done
echo "----"; echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]

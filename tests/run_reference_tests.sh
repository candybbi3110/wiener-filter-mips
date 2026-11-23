#!/usr/bin/env bash
# Runs the C++ reference model on every test case and compares with expected_output.txt.
# Usage: tests/run_reference_tests.sh [path/to/wiener_ref]   (default: build/wiener_ref)
set -u
BIN="${1:-build/wiener_ref}"
DIR="$(cd "$(dirname "$0")" && pwd)"
pass=0; fail=0
for tc in "$DIR"/tc*/; do
  name="$(basename "$tc")"
  out="$(mktemp)"
  "$BIN" "$tc/input.txt" "$tc/desired.txt" "$out" > /dev/null
  # whitespace-insensitive compare (MARS output has trailing spaces)
  if diff -wB "$out" "$tc/expected_output.txt" > /dev/null; then
    echo "PASS  $name"; pass=$((pass+1))
  else
    echo "FAIL  $name"; diff -wB "$out" "$tc/expected_output.txt"; fail=$((fail+1))
  fi
  rm -f "$out"
done
echo "----"; echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]

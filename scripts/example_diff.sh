#!/usr/bin/env bash
# Equivalence harness for backend and renderer changes.
#
# The test suite asserts on plotter and renderer call sequences, not on
# serialised backend output, so it misses output regressions (plotly trace
# order, bokeh array serialisation, matplotlib spines). This script
# regenerates examples/ and lists the files whose content differs from HEAD
# after normalising generated element ids. Matplotlib PNGs are compared
# byte for byte.
#
# Usage: scripts/example_diff.sh [--keep]
#   --keep  leave the regenerated files in the working tree (to commit an
#           intended change); by default they are restored.
# Prints NO REAL DIFFS or one REAL DIFF: <file> line per changed export.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
keep="${1:-}"
uv run python examples/generate_example_plots.py >/tmp/example-gen.log 2>&1 || {
  echo "GENERATOR FAILED"; tail -30 /tmp/example-gen.log; exit 2; }
norm() {
  perl -pe '
    s/[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/UUID/g;
    s/\b[0-9a-f]{32}\b/HEX32/g;
    s/\bp[0-9]{3,6}\b/PID/g;
  '
}
real=()
while IFS= read -r f; do
  case "$f" in
    *.html)
      if ! diff -q <(git show "HEAD:$f" | norm) <(norm < "$f") >/dev/null; then real+=("$f"); fi ;;
    *) real+=("$f") ;;
  esac
done < <(git diff --name-only -- examples/)
if [ "${#real[@]}" -eq 0 ]; then echo "NO REAL DIFFS"; else printf 'REAL DIFF: %s\n' "${real[@]}"; fi
if [ "$keep" != "--keep" ] && [ -n "$(git diff --name-only -- examples/)" ]; then
  git stash push -- examples/ >/dev/null && git stash drop >/dev/null
fi

#!/usr/bin/env bash
# Run @probelabs/maid over markdown files and fail loudly on real errors.
#
# We treat maid's exit code as the source of truth:
#   0 = clean
#   1 = lint errors found (we still inspect output to filter
#       FL-STYLE-TARGET-UNKNOWN false positives)
#   anything else = maid crashed (npx/network/OOM/module-not-found);
#                   surface the output and abort.
#
# grep -c returns exit 1 on zero matches; under `set -e` that aborts the
# script, so we capture with `|| real_errors=0` and then assert the value
# is numeric before comparison.
set -euo pipefail

fail=0
for f in "$@"; do
  set +e
  output=$(npx -y @probelabs/maid@0.0.29 "$f" 2>&1)
  rc=$?
  set -e
  if [ "$rc" -ne 0 ] && [ "$rc" -ne 1 ]; then
    echo "ERROR: maid crashed (exit $rc) on $f" >&2
    echo "$output" >&2
    exit 2
  fi
  # FL-STYLE-TARGET-UNKNOWN: maid false-positive on some class-diagram
  # style references; safe to filter until upstream fixes it.
  real_errors=$(printf '%s\n' "$output" | grep -v 'FL-STYLE-TARGET-UNKNOWN' | grep -cE '\berror\b') \
    || real_errors=0
  case "$real_errors" in
    ''|*[!0-9]*)
      echo "ERROR: mermaid-lint produced non-numeric error count on $f" >&2
      echo "$output" >&2
      exit 2
      ;;
  esac
  if [ "$real_errors" -gt 0 ]; then
    echo "$f:"
    echo "$output"
    fail=1
  fi
done
exit "$fail"

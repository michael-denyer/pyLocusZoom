#!/usr/bin/env bash
# Run @probelabs/maid over markdown files and fail loudly on real errors.
#
# grep -c returns exit 1 when there are zero matches; under `set -e` that
# aborts the script. Capture the count with `|| real_errors=0` so the
# zero-match case is distinguishable from a grep crash or broken pipe.
# Then assert the value is numeric before comparison so a silent empty
# string can't sneak through.
set -euo pipefail

fail=0
for f in "$@"; do
  output=$(npx -y @probelabs/maid@0.0.29 "$f" 2>&1) || true
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

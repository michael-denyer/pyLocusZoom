#!/usr/bin/env bash
# Check isolated exports against HEAD. --keep accepts changes to clean exports.
set -euo pipefail
exec python3 "$(dirname "$0")/example_diff.py" "$@"

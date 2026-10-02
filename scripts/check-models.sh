#!/usr/bin/env bash
# Check the formal models: every specs/tla/*.matrix through TLC and the
# specs/lean package through the Lean checker. The runners are the helper
# scripts of https://github.com/michael-denyer/agent-formal-verify, so
# FORMAL_VERIFY must name a checkout of it; CI pins the commit in ci.yml.
#
# Usage: FORMAL_VERIFY=<checkout> scripts/check-models.sh [tla|lean]
#   JAVA=<path>  the java binary TLC runs on (default: java on PATH)
set -euo pipefail
usage="usage: FORMAL_VERIFY=<agent-formal-verify checkout> $0 [tla|lean]"
helpers=${FORMAL_VERIFY:?$usage}/skills/formal-verify/scripts
[ -f "$helpers/setup.sh" ] || { echo "FAIL $helpers/setup.sh not found; $usage"; exit 1; }
root=$(cd "$(dirname "$0")/.." && pwd)

check_tla() {
  local status=0 matrix
  bash "$helpers/setup.sh" tla || return 1
  for matrix in "$root"/specs/tla/*.matrix; do
    bash "$helpers/tlc-matrix.sh" "$matrix" || status=1
  done
  return "$status"
}

check_lean() {
  bash "$helpers/setup.sh" lean "$root/specs/lean" || return 1
  bash "$helpers/lean-check.sh" "$root/specs/lean"
}

case ${1:-all} in
  tla) check_tla ;;
  lean) check_lean ;;
  all)
    status=0
    check_tla || status=1
    check_lean || status=1
    exit "$status"
    ;;
  *) echo "$usage"; exit 2 ;;
esac

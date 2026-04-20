#!/usr/bin/env bash
# Render-parity check: run mermaid-cli (the real renderer) over each
# markdown file to catch constructs maid accepts but mermaid-js rejects.
#
# Puppeteer's bundled Chromium can't use the user-namespace sandbox on
# GitHub's Ubuntu runners (AppArmor restriction) and on many local Linux
# machines. Configure mmdc to launch Chromium with --no-sandbox.
set -euo pipefail

# Use a scratch dir so we control extensions (mmdc requires .svg/.png/.pdf
# suffix) and avoid mktemp -t portability differences between macOS and Linux.
tmpdir=$(mktemp -d -t mmdc-check.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
config_file="$tmpdir/puppeteer-config.json"
out_svg="$tmpdir/out.svg"
printf '{"args": ["--no-sandbox"]}\n' > "$config_file"

fail=0
for f in "$@"; do
  # mmdc writes one SVG per diagram as out-1.svg, out-2.svg, ... next to
  # $out_svg. Clear prior run's outputs so we can validate this file's.
  rm -f "$tmpdir"/out*.svg
  set +e
  mmdc_output=$(npx -y @mermaid-js/mermaid-cli@11.4.2 \
      -p "$config_file" \
      -i "$f" -o "$out_svg" 2>&1)
  rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then
    echo "ERROR: mmdc failed on $f (exit $rc)" >&2
    echo "$mmdc_output" >&2
    fail=1
    continue
  fi
  # mmdc can exit 0 while producing an empty or invalid SVG under broken
  # puppeteer/Chromium configs. Verify the outputs actually rendered.
  # For .md inputs mmdc writes out-1.svg, out-2.svg, ...; for non-markdown
  # it writes $out_svg directly.
  shopt -s nullglob
  produced=("$tmpdir"/out*.svg)
  shopt -u nullglob
  # Files with no mermaid blocks legitimately produce no output.
  # mmdc prints "No mermaid charts found" in that case — treat as success.
  if [ "${#produced[@]}" -eq 0 ]; then
    if printf '%s' "$mmdc_output" | grep -q 'No mermaid charts found'; then
      continue
    fi
    echo "ERROR: mmdc produced no SVG output for $f" >&2
    echo "$mmdc_output" >&2
    fail=1
    continue
  fi
  for svg in "${produced[@]}"; do
    if [ ! -s "$svg" ] || ! grep -q '<svg' "$svg"; then
      echo "ERROR: mmdc produced empty or invalid SVG ($svg) for $f" >&2
      echo "$mmdc_output" >&2
      fail=1
      break
    fi
  done
done
exit "$fail"

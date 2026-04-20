#!/usr/bin/env bash
# Render-parity check: run mermaid-cli (the real renderer) over each
# markdown file to catch constructs maid accepts but mermaid-js rejects.
#
# Puppeteer's bundled Chromium can't use the user-namespace sandbox on
# GitHub's Ubuntu runners (AppArmor restriction) and on many local Linux
# machines. Configure mmdc to launch Chromium with --no-sandbox.
set -euo pipefail

config_file=$(mktemp -t puppeteer-config.XXXXXX.json)
trap 'rm -f "$config_file"' EXIT
printf '{"args": ["--no-sandbox"]}\n' > "$config_file"

fail=0
for f in "$@"; do
  if ! npx -y @mermaid-js/mermaid-cli@11.4.2 \
      -p "$config_file" \
      -i "$f" -o /tmp/mmdc_check.svg; then
    echo "mmdc failed on $f" >&2
    fail=1
  fi
done
exit "$fail"

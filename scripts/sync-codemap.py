"""Keep docs/CODEMAP.md anchor line numbers in sync with source markers.

Every row in the CODEMAP tables is expected to be backed by a per-symbol anchor
comment of the form ``# [<id>:<symbol>]`` on the line immediately above the
target ``def``/``class``. The script does two things:

1. Verifies that each CODEMAP row's claimed line is the line *below* an anchor
   comment for the same ``<id>:<symbol>`` pair.
2. In ``--fix`` mode, rewrites the table's ``file:line`` reference to match the
   anchor's actual location in the source tree.

Exit status:
- ``0`` — CODEMAP is in sync (check mode) or was updated successfully
  (``--fix`` mode).
- ``1`` — drift detected in check mode, or an anchor is missing entirely.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CODEMAP = ROOT / "docs" / "CODEMAP.md"
SRC_ROOT = ROOT / "src" / "pylocuszoom"

# Matches `| 1a | LocusZoomPlotter | ... | [plotter.py:58](../src/pylocuszoom/plotter.py#L58) |`
ROW_RE = re.compile(
    r"^(?P<prefix>\|\s*(?P<id>\d[a-z])\s*\|\s*(?P<symbol>[A-Za-z_][\w]*)\s*\|[^|]*\|\s*)"
    r"\[(?P<label>[^\]]+)\]\((?P<href>[^)]+)\)"
    r"(?P<suffix>\s*\|\s*)$"
)

ANCHOR_RE = re.compile(r"#\s*\[(?P<id>\d[a-z]):(?P<symbol>[A-Za-z_][\w]*)\]")


def find_anchor(id_: str, symbol: str) -> tuple[Path, int] | None:
    """Locate the line number of the symbol referenced by ``# [id:symbol]``."""
    token = f"[{id_}:{symbol}]"
    for path in SRC_ROOT.rglob("*.py"):
        with path.open() as fh:
            lines = fh.readlines()
        for idx, line in enumerate(lines):
            if token not in line:
                continue
            # Symbol line is the first non-blank, non-comment line after the anchor.
            for follow_idx in range(idx + 1, len(lines)):
                stripped = lines[follow_idx].strip()
                if not stripped or stripped.startswith("#"):
                    continue
                return path, follow_idx + 1  # 1-indexed
            return path, idx + 2
    return None


def process(fix: bool) -> int:
    text = CODEMAP.read_text()
    new_lines: list[str] = []
    errors: list[str] = []
    updates = 0

    for lineno, line in enumerate(text.splitlines(keepends=True), start=1):
        match = ROW_RE.match(line.rstrip("\n"))
        if not match:
            new_lines.append(line)
            continue

        id_ = match["id"]
        symbol = match["symbol"]
        found = find_anchor(id_, symbol)
        if found is None:
            errors.append(
                f"CODEMAP.md:{lineno} references [{id_}:{symbol}] but no anchor "
                f"comment '# [{id_}:{symbol}]' exists under src/pylocuszoom/."
            )
            new_lines.append(line)
            continue

        anchor_path, anchor_line = found
        rel_path = anchor_path.relative_to(ROOT).as_posix()
        label = f"{anchor_path.name}:{anchor_line}"
        href = f"../{rel_path}#L{anchor_line}"
        expected_link = f"[{label}]({href})"

        current_link = f"[{match['label']}]({match['href']})"
        if current_link == expected_link:
            new_lines.append(line)
            continue

        updated = f"{match['prefix']}{expected_link}{match['suffix']}\n"
        new_lines.append(updated)
        if fix:
            updates += 1
        else:
            errors.append(
                f"CODEMAP.md:{lineno} points at {current_link} "
                f"but anchor is at {expected_link}"
            )

    if errors:
        for err in errors:
            print(err, file=sys.stderr)
        if not fix:
            print(
                f"\n{len(errors)} drift(s) detected. Run `python scripts/sync-codemap.py --fix`.",
                file=sys.stderr,
            )
            return 1

    if fix and updates:
        CODEMAP.write_text("".join(new_lines))
        print(f"Updated {updates} CODEMAP.md row(s).")
    elif fix:
        print("CODEMAP.md already in sync.")

    return 1 if errors and not fix else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fix", action="store_true", help="Rewrite drifted rows.")
    args = parser.parse_args()
    return process(fix=args.fix)


if __name__ == "__main__":
    sys.exit(main())

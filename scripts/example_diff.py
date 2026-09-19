#!/usr/bin/env python3
"""Compare isolated example exports with HEAD; explicitly accept with --keep."""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

EXPORT_TYPES = {".html", ".png"}


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args])


def normalized(data, path):
    if path.suffix != ".html":
        return data
    for pattern, replacement in (
        (rb"[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}", b"UUID"),
        (rb"\b[0-9a-f]{32}\b", b"HEX32"),
        (rb"\bp[0-9]{3,6}\b", b"PID"),
    ):
        data = re.sub(pattern, replacement, data)
    return data


def compare(repo, generated):
    tracked = {
        Path(os.fsdecode(name))
        for name in git(
            repo, "ls-tree", "-rz", "--name-only", "HEAD", "--", "examples/"
        ).split(b"\0")
        if name and Path(os.fsdecode(name)).suffix in EXPORT_TYPES
    }
    candidates = {
        p.relative_to(generated)
        for p in (generated / "examples").rglob("*")
        if p.is_file() and p.suffix in EXPORT_TYPES
    }
    changed = []
    for path in sorted(tracked | candidates):
        before = (
            git(repo, "show", f"HEAD:{path.as_posix()}") if path in tracked else None
        )
        after = (generated / path).read_bytes() if path in candidates else None
        if (
            before is None
            or after is None
            or normalized(before, path) != normalized(after, path)
        ):
            changed.append(path)
    return changed, tracked


def accept(repo, generated, changed, tracked):
    # Check every destination before touching any of them. Manual export edits
    # must survive even when the caller deliberately accepts generated output.
    dirty = {
        Path(os.fsdecode(p))
        for p in git(
            repo, "diff", "HEAD", "--name-only", "-z", "--", "examples/"
        ).split(b"\0")
        if p
    }
    conflicts = [
        p for p in changed if p in dirty or (p not in tracked and (repo / p).exists())
    ]
    if conflicts:
        raise ValueError(
            "Refusing to overwrite modified exports: " + ", ".join(map(str, conflicts))
        )
    for path in changed:
        source, target = generated / path, repo / path
        if source.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        else:
            target.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep",
        action="store_true",
        help="accept changed exports, refusing manually modified files",
    )
    options = parser.parse_args()
    repo = Path(
        os.fsdecode(
            subprocess.check_output(["git", "rev-parse", "--show-toplevel"])
        ).strip()
    )
    with tempfile.TemporaryDirectory(prefix="pylocuszoom-examples-") as directory:
        generated = Path(directory)
        for backend in ("matplotlib", "plotly", "bokeh"):
            (generated / "examples" / backend).mkdir(parents=True)
        environment = dict(os.environ)
        environment["PYTHONPATH"] = (
            str(repo / "src") + os.pathsep + environment.get("PYTHONPATH", "")
        )
        command = [
            "uv",
            "run",
            "--project",
            str(repo),
            "python",
            str(repo / "examples/generate_example_plots.py"),
        ]
        with (generated / "generator.log").open("w+") as log:
            result = subprocess.run(
                command,
                cwd=generated,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if result.returncode:
                log.seek(0)
                print("GENERATOR FAILED", file=sys.stderr)
                print("".join(log.readlines()[-30:]), file=sys.stderr)
                return 2
        changed, tracked = compare(repo, generated)
        for path in changed:
            print(f"REAL DIFF: {path}")
        if not changed:
            print("NO REAL DIFFS")
        if options.keep:
            accept(repo, generated, changed, tracked)
            return 0
        return int(bool(changed))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(error, file=sys.stderr)
        sys.exit(2)

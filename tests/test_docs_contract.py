"""Guard the two documentation facts that drifted from the code.

Both checks exist because prose cannot be trusted to track a value that lives
somewhere else. The first pins every pytest command line in the docs to what
``addopts`` already supplies. The second pins the public-surface table in
CODEMAP.md to ``pylocuszoom.__all__``.
"""

import re
from pathlib import Path

import tomllib

import pylocuszoom

REPO_ROOT = Path(__file__).resolve().parents[1]

DOCUMENTED_COMMAND_FILES = [
    "CONTRIBUTING.md",
    "README.md",
    "docs/TESTING.md",
    "docs/DEVELOPMENT.md",
    ".pre-commit-config.yaml",
    ".github/workflows/ci.yml",
]

INVOCATION = re.compile(r"(?:uv run|run:|entry:)[^|]*pytest")

REDUNDANT_FLAGS = [
    re.compile(r"-n\s+\d"),
    re.compile(r"--timeout[=\s]"),
    re.compile(r"--cov(?![a-z-])"),
    re.compile(r"--cov-report"),
    re.compile(r"(?<![\w-])-v(?![\w-])"),
    re.compile(r'-m\s+["\']?not\s+integration'),
]


def _addopts() -> str:
    config = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    return config["tool"]["pytest"]["ini_options"]["addopts"]


def test_addopts_still_supplies_every_flag_this_guard_strips():
    """The flags this test bans are exactly the ones addopts sets."""
    addopts = _addopts()

    for flag in (
        "-n 3",
        "--timeout=",
        "--cov=",
        "--cov-report=",
        "-v",
        "not integration",
    ):
        assert flag in addopts


def test_no_doc_repeats_a_flag_addopts_already_supplies():
    """A pytest command in the docs carries no flag addopts already sets."""
    offenders = []

    for relative in DOCUMENTED_COMMAND_FILES:
        path = REPO_ROOT / relative
        for number, line in enumerate(path.read_text().splitlines(), start=1):
            if not INVOCATION.search(line) or "--no-cov" in line:
                continue
            for pattern in REDUNDANT_FLAGS:
                if pattern.search(line):
                    offenders.append(f"{relative}:{number}: {line.strip()}")
                    break

    assert not offenders, (
        "These command lines repeat a flag that addopts already supplies:\n"
        + "\n".join(offenders)
    )


def test_codemap_lists_every_public_name():
    """Every name in __all__ appears in the CODEMAP public-surface table."""
    codemap = (REPO_ROOT / "docs" / "CODEMAP.md").read_text()
    missing = [name for name in pylocuszoom.__all__ if f"`{name}`" not in codemap]

    assert not missing, (
        "Add these to the Public API Surface table in docs/CODEMAP.md: "
        + ", ".join(sorted(missing))
    )

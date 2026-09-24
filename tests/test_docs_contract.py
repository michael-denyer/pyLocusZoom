"""Guard the documentation facts that drifted from the code.

Each check exists because prose cannot be trusted to track a value that lives
somewhere else. The first pins every pytest command line in the docs to what
``addopts`` already supplies. The second pins the public-surface table in
CODEMAP.md to ``pylocuszoom.__all__``. The third pins the USER_GUIDE API
Stability tables to the core and toolbox tiers of ``__all__``.
"""

import re
from pathlib import Path

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


def test_addopts_still_supplies_every_flag_this_guard_strips(pytestconfig):
    """The flags this test bans are exactly the ones addopts sets."""
    addopts = " ".join(pytestconfig.getini("addopts"))

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


def _all_tiers() -> tuple[set[str], set[str]]:
    """Split ``__all__`` in ``__init__.py`` at its ``# Toolbox`` comment."""
    source = (REPO_ROOT / "src" / "pylocuszoom" / "__init__.py").read_text()
    block = source[source.index("__all__ = [") :]
    block = block[: block.index("\n]")]
    core, toolbox = block.split("# Toolbox", 1)
    return set(re.findall(r'"(\w+)"', core)), set(re.findall(r'"(\w+)"', toolbox))


def _guide_tiers() -> tuple[set[str], set[str]]:
    """Read the names in the USER_GUIDE API Stability Core and Toolbox tables."""
    guide = (REPO_ROOT / "docs" / "USER_GUIDE.md").read_text()
    section = guide[guide.index("## API Stability") :]
    section = section[: section.index("\n---")]
    core, toolbox = section.split("### Toolbox", 1)

    def table_names(text: str) -> set[str]:
        rows = (line for line in text.splitlines() if line.startswith("|"))
        return {name for row in rows for name in re.findall(r"`(\w+)`", row)}

    return table_names(core), table_names(toolbox)


def test_all_tiers_partition_the_exported_names():
    """The two commented blocks of __all__ hold every export exactly once."""
    core, toolbox = _all_tiers()

    assert not core & toolbox
    assert core | toolbox == set(pylocuszoom.__all__)


def test_user_guide_api_stability_tables_match_all_tiers():
    """The USER_GUIDE Core and Toolbox tables list exactly the two __all__ tiers."""
    core, toolbox = _all_tiers()
    guide_core, guide_toolbox = _guide_tiers()

    assert guide_core == core, (
        f"Core tier missing from USER_GUIDE: {sorted(core - guide_core)}; "
        f"listed but not core: {sorted(guide_core - core)}"
    )
    assert guide_toolbox == toolbox, (
        f"Toolbox tier missing from USER_GUIDE: {sorted(toolbox - guide_toolbox)}; "
        f"listed but not toolbox: {sorted(guide_toolbox - toolbox)}"
    )

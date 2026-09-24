"""Example verification never uses the checkout as a generation workspace."""

import os
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/example_diff.sh"


@pytest.fixture
def example_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    # Git hooks export repository-specific state; this fixture owns another repo.
    environment = dict(
        {key: value for key, value in os.environ.items() if not key.startswith("GIT_")},
        GIT_AUTHOR_NAME="Test",
        GIT_AUTHOR_EMAIL="test@example.invalid",
        GIT_COMMITTER_NAME="Test",
        GIT_COMMITTER_EMAIL="test@example.invalid",
    )

    def git(*args):
        return subprocess.run(
            ["git", *args],
            cwd=repo,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        ).stdout

    git("init", "-q")
    (repo / "examples/matplotlib").mkdir(parents=True)
    (repo / "examples/matplotlib/plot.png").write_bytes(b"baseline")
    (repo / "examples/generate_example_plots.py").write_text("# baseline generator\n")
    git("add", ".")
    git(
        "-c",
        "commit.gpgsign=false",
        "-c",
        "core.hooksPath=/dev/null",
        "commit",
        "-qm",
        "baseline",
    )
    binaries = tmp_path / "bin"
    binaries.mkdir()
    uv = binaries / "uv"
    uv.write_text("#!/bin/sh\nprintf generated > examples/matplotlib/plot.png\n")
    uv.chmod(0o700)
    environment["PATH"] = str(binaries) + os.pathsep + environment["PATH"]
    return repo, environment, git


def run_check(fixture, *args):
    repo, environment, _ = fixture
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        cwd=repo,
        env=environment,
        capture_output=True,
        text=True,
    )


def test_check_preserves_manual_edits_and_returns_difference(example_repo):
    repo, _, git = example_repo
    generator = repo / "examples/generate_example_plots.py"
    generator.write_text("# manual work\n")
    result = run_check(example_repo)
    assert generator.read_text() == "# manual work\n"
    assert (repo / "examples/matplotlib/plot.png").read_bytes() == b"baseline"
    assert result.returncode == 1
    assert result.stdout.strip() == "REAL DIFF: examples/matplotlib/plot.png"
    assert git("stash", "list") == ""


def test_accept_updates_only_generated_exports(example_repo):
    repo, _, _ = example_repo
    generator = repo / "examples/generate_example_plots.py"
    generator.write_text("# manual work\n")
    result = run_check(example_repo, "--keep")
    assert result.returncode == 0, result.stderr
    assert generator.read_text() == "# manual work\n"
    assert (repo / "examples/matplotlib/plot.png").read_bytes() == b"generated"


def test_accept_refuses_to_overwrite_modified_export(example_repo):
    repo, _, _ = example_repo
    plot = repo / "examples/matplotlib/plot.png"
    plot.write_bytes(b"manual export")
    result = run_check(example_repo, "--keep")
    assert result.returncode == 2
    assert plot.read_bytes() == b"manual export"


def test_identical_exports_return_success(example_repo):
    _, environment, _ = example_repo
    uv = Path(environment["PATH"].split(os.pathsep)[0]) / "uv"
    uv.write_text("#!/bin/sh\nprintf baseline > examples/matplotlib/plot.png\n")
    result = run_check(example_repo)
    assert result.returncode == 0
    assert result.stdout.strip() == "NO REAL DIFFS"


@pytest.mark.parametrize("args", [(), ("--keep",)])
def test_generator_failure_preserves_checkout(example_repo, args):
    # A generator that fails after writing some exports (as it does on a
    # degraded figure) must leave the checkout untouched, even with --keep.
    repo, environment, git = example_repo
    uv = Path(environment["PATH"].split(os.pathsep)[0]) / "uv"
    uv.write_text(
        "#!/bin/sh\nprintf incomplete > examples/matplotlib/plot.png\nexit 3\n"
    )
    result = run_check(example_repo, *args)
    assert result.returncode == 2
    assert "GENERATOR FAILED" in result.stderr
    assert (repo / "examples/matplotlib/plot.png").read_bytes() == b"baseline"
    assert git("status", "--porcelain") == ""


CDN_HTML = (
    '<script charset="utf-8" src="https://cdn.plot.ly/plotly-{version}.min.js">'
    '</script><div id="{uuid}"></div>'
)


@pytest.mark.parametrize(
    ("generated_version", "expected"),
    [
        ("3.5.0", "NO REAL DIFFS"),
        ("3.6.0", "REAL DIFF: examples/plotly/plot.html"),
    ],
)
def test_cdn_plotly_export_compares_by_plotly_version(
    example_repo, generated_version, expected
):
    # A regeneration changes the div UUID but not the CDN tag; a plotly.js
    # upgrade changes the runtime the export loads, so it is a real change.
    repo, environment, git = example_repo
    (repo / "examples/plotly").mkdir()
    (repo / "examples/plotly/plot.html").write_text(
        CDN_HTML.format(version="3.5.0", uuid="00000000-0000-0000-0000-000000000000")
        + "\n"
    )
    git("add", ".")
    git("-c", "core.hooksPath=/dev/null", "commit", "-qm", "html baseline")
    generated = CDN_HTML.format(
        version=generated_version, uuid="11111111-1111-1111-1111-111111111111"
    )
    uv = Path(environment["PATH"].split(os.pathsep)[0]) / "uv"
    uv.write_text(
        "#!/bin/sh\nprintf baseline > examples/matplotlib/plot.png\n"
        f"cat > examples/plotly/plot.html <<'HTML'\n{generated}\nHTML\n"
    )
    result = run_check(example_repo)
    assert result.stdout.strip() == expected

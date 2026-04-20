# Contributing to pyLocusZoom

Thanks for your interest in contributing. This guide covers how to set up a
development environment, the standards your changes must meet, and how to get
a pull request merged.

## Development setup

See the [README](README.md) for installation instructions and
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for a tour of the codebase. For
configuration knobs (species flags, recombination caches, LD reference files),
see [docs/CONFIGURATION.md](docs/CONFIGURATION.md).

Short version for contributors:

```bash
git clone https://github.com/michael-denyer/pyLocusZoom.git
cd pyLocusZoom
uv sync
uv tool install prek && prek install   # Rust pre-commit, ~10x faster than pre-commit
uv run python -m pytest tests/ -n 3
```

The hook install step is required — without it, you will only discover lint,
format, markdown, or mermaid failures in CI. We recommend [prek](https://github.com/j178/prek),
a drop-in replacement for `pre-commit` written in Rust; it reads the same
`.pre-commit-config.yaml` but installs and runs hooks dramatically faster. If
you already use `pre-commit`, `uv run pre-commit install` still works. The
mermaid hooks shell out to `npx`, so you need Node.js 20+ on your PATH
(`brew install node` on macOS). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)
for the full dev setup.

Python 3.10, 3.11, and 3.12 are all supported and tested in CI.

## Coding standards

- **Formatter and linter:** [ruff](https://docs.astral.sh/ruff/), pinned to
  version `0.15.2` in CI. Run both before committing:
  ```bash
  uv tool run ruff check src/ tests/
  uv tool run ruff format --check src/ tests/
  ```
  If the format check fails, run `uv tool run ruff format src/ tests/` to fix.
- **Tests:** `pytest` with `pytest-xdist` for parallelism, `pytest-randomly`
  for order-independence, and `pytest-timeout` for hang detection. Run:
  ```bash
  uv run python -m pytest tests/ -n 3
  ```
- **Docstrings:** Google-style for all public functions and classes.
- **Mermaid diagrams:** validated in CI with both `@probelabs/maid` (syntax)
  and `@mermaid-js/mermaid-cli` (renderer parity). Keep diagrams in fenced
  ` ```mermaid ` blocks.
- **Markdown:** `markdownlint-cli2` runs on every `.md` file in CI.
- **Links:** `lychee` checks every link in committed markdown; broken links
  fail the build.

CI enforces all of the above on every pull request via
[`.github/workflows/ci.yml`](.github/workflows/ci.yml). A PR with lint, format,
markdown, mermaid, link-check, or test failures will not merge.

## Pull request guidelines

- **Base branch:** open PRs against `main`.
- **Tests first:** follow test-driven development — add or update tests in
  `tests/` before (or alongside) the implementation. Mock PLINK calls rather
  than requiring a local install; see `tests/test_ld.py` for the pattern.
- **Changelog:** add an entry to `CHANGELOG.md` under the `## [Unreleased]`
  section, using the `Added` / `Changed` / `Fixed` / `Removed` categories.
- **Docs:** update `README.md`, `docs/USER_GUIDE.md`, `docs/ARCHITECTURE.md`,
  or `docs/CONFIGURATION.md` when your change alters public APIs, adds
  features, or changes behavior that users rely on.
- **Example plots:** if your change affects rendering, regenerate them with
  `uv run python examples/generate_example_plots.py` and commit the updated
  images.
- **Commits:** keep messages focused on *what* changed and *why*. Do not
  include AI or tool attribution.
- **Scope:** one logical change per PR. Refactors and feature work belong in
  separate PRs.

## Issue reporting

Bugs and feature requests go through GitHub Issues:
<https://github.com/michael-denyer/pyLocusZoom/issues>

Templates are provided and should be used:

- **Bug reports**
  ([`.github/ISSUE_TEMPLATE/bug_report.md`](.github/ISSUE_TEMPLATE/bug_report.md))
  must include a minimal reproducible example, expected vs actual behavior,
  and your environment (pylocuszoom version, Python version, OS, install
  method).
- **Feature requests**
  ([`.github/ISSUE_TEMPLATE/feature_request.md`](.github/ISSUE_TEMPLATE/feature_request.md))
  should describe the use case and proposed API, not just the feature name.

Before opening an issue, check [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md) —
many common questions are answered there.

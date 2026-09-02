# Development

Guidance for contributors working on pyLocusZoom locally. For first-run install and
quickstart, see [GETTING-STARTED.md](GETTING-STARTED.md) or the
[README](../README.md). For system architecture, see [ARCHITECTURE.md](ARCHITECTURE.md).
For environment variables and configuration, see [CONFIGURATION.md](CONFIGURATION.md).

## Local Setup

pyLocusZoom uses [`uv`](https://docs.astral.sh/uv/) for dependency management and
virtualenv handling. Python 3.10, 3.11, or 3.12 is required (`requires-python = ">=3.10"`
in `pyproject.toml`).

1. Fork the repository on GitHub, then clone your fork:

   ```bash
   git clone https://github.com/<your-username>/pyLocusZoom.git
   cd pyLocusZoom
   ```

2. Install `uv` if you do not have it: see [uv installation docs](https://docs.astral.sh/uv/getting-started/installation/).

3. Create the virtualenv and install all dependencies (runtime, dev, and optional extras):

   ```bash
   uv sync --all-extras
   ```

   This reads `pyproject.toml` and `uv.lock` to produce a reproducible `.venv/`.

4. Install Node.js 20+ on your PATH. The `mermaid-lint` and `mermaid-render` pre-commit
   hooks call `npx`, and CI runs them on every PR. Without Node, pre-commit will fail
   with `npx: command not found`.
   - macOS: `brew install node`
   - Ubuntu/Debian: `sudo apt install nodejs npm` (or use [nvm](https://github.com/nvm-sh/nvm))

5. Install the pre-commit hooks so lint/format/tests run automatically on `git commit`:

   ```bash
   uv run pre-commit install
   ```

6. Verify the setup:

   ```bash
   uv run python -m pytest tests/ -n 3
   ```

## Build Commands

pyLocusZoom has no `npm`-style `scripts` block; commands are run directly via `uv`.
The common development commands are:

| Command | Description |
|---------|-------------|
| `uv sync --all-extras` | Install all dependencies (runtime + `dev` + `spark` + `all` extras) into `.venv/`. |
| `uv run python -m pytest tests/ -n 3` | Run the full test suite in parallel across 3 workers (matches CI). |
| `uv run python -m pytest tests/ -n 3 --no-cov` | Fast iteration: skip coverage reporting. |
| `uv run python -m pytest tests/test_plotter.py -v` | Run a single test file. |
| `uv run python -m pytest tests/ -m "not integration"` | Skip integration tests (default in CI and `pyproject.toml` addopts). |
| `uv tool run ruff check src/ tests/` | Run ruff lint checks (no fixes). |
| `uv tool run ruff format src/ tests/` | Apply ruff formatting. |
| `uv tool run ruff format --check src/ tests/` | Verify formatting without changes (matches CI check). |
| `uv run pre-commit run --all-files` | Run the full pre-commit suite against every file in the repo. |
| `uv build` | Build the wheel and sdist via hatchling into `dist/`. |
| `uv run python examples/generate_example_plots.py` | Regenerate example plots shown in the README. |
| `scripts/example_diff.sh [--keep]` | Regenerate the examples and list the exports whose content changed after id normalisation; the equivalence check for backend and renderer changes. |
| `uv lock` | Refresh `uv.lock` after changing dependencies in `pyproject.toml`. |

See [CONTRIBUTING.md](../CONTRIBUTING.md) for the full pre-commit and pre-PR checklists.

## Code Style

Formatting and linting are handled by a **single tool: [ruff](https://github.com/astral-sh/ruff)**.
Configuration lives in `pyproject.toml` under `[tool.ruff]` and `[tool.ruff.lint]`.

- Line length: **88 characters** (`line-length = 88`)
- Target Python: `py310` (`target-version = "py310"`)
- Enabled lint rule sets: `E`, `F`, `I`, `W` (pycodestyle errors/warnings, pyflakes, isort)
- Ignored rules: `E501` (line-length handled by the formatter, not the linter)
- Per-file overrides: `examples/*` allows `E402` (mid-file imports, useful for tutorial-style scripts)

Run the checks locally before pushing:

```bash
uv tool run ruff check src/ tests/
uv tool run ruff format --check src/ tests/
```

The ruff version is pinned to `0.15.2` in CI (`.github/workflows/ci.yml`) and in the
pre-commit config (`.pre-commit-config.yaml`). Keep local `uv tool` usage in sync to
avoid "works locally, fails in CI" drift.

### Docstrings

Use **Google-style docstrings** for all public functions and classes. Key conventions:

- First line is a short imperative summary.
- Do not duplicate type annotations in the docstring — they belong on the signature.
- Include an `Example` section for non-trivial public APIs.
- Omit sections (`Raises`, `Example`) that do not apply.

For a concrete example, see `prepare_manhattan_data` in
[`src/pylocuszoom/manhattan.py`](../src/pylocuszoom/manhattan.py).

### Pre-commit Hooks

`.pre-commit-config.yaml` wires up every quality check that CI runs. The active hooks are:

| Hook | Purpose |
|------|---------|
| `pre-commit-hooks` v5.0.0 (check-yaml, check-toml, check-json, check-ast, check-merge-conflict, check-added-large-files, end-of-file-fixer, trailing-whitespace) | Generic hygiene checks |
| `ruff` v0.15.2 (`--fix --exit-non-zero-on-fix`) | Lint and auto-fix |
| `ruff-format` v0.15.2 | Format |
| `lychee` v0.23.0 | Markdown link checking (`lychee.toml`) |
| `markdownlint-cli2` 0.22.0 | Markdown style |
| `mermaid-lint` (`scripts/check-mermaid-maid.sh`) | Fast mermaid syntax validation via `@probelabs/maid` |
| `mermaid-render` (`scripts/check-mermaid-render.sh`) | Renderer-parity validation via `mmdc` (catches constructs maid accepts but the real renderer rejects) |
| `yamllint` 1.37.0 | YAML linting (`.yamllint.yaml`) |
| `no-planning-files` | Blocks accidental commits of `.planning/` |
| `no-gitignored-files` | Blocks staging of gitignored files |
| `pytest-cov` | Runs `uv run python -m pytest -n 3` on every Python change |

## Branch Conventions

- The default branch is **`main`**. CI runs on `push` to `main` and on all pull requests
  targeting `main` (`.github/workflows/ci.yml`).
- Feature branches follow the pattern shown in [CONTRIBUTING.md](../CONTRIBUTING.md#pull-request-process):
  `git checkout -b feature/your-feature`. Use `fix/…` for bug fixes and `chore/…` for
  maintenance tasks. No stricter convention is enforced.
- **Always commit and push before ending a session.** Do not leave `uv.lock` or other
  regenerated files uncommitted after a dependency bump.

## PR Process

1. Work on a feature branch off `main` (see above).
2. Write tests **before** the implementation where practical
   (see [docs/TESTING.md](TESTING.md) on test-driven development).
3. Before opening the PR, run the pre-commit checklist locally:

   ```bash
   uv run python -m pytest tests/ -n 3
   uv tool run ruff check src/ tests/
   uv tool run ruff format --check src/ tests/
   ```

4. Update [CHANGELOG.md](../CHANGELOG.md) under the `## [Unreleased]` section. Every
   PR is expected to add an entry under `Added`, `Changed`, `Fixed`, or `Removed`.
5. If your change alters architecture, data flow, or adds a new module, update the
   mermaid diagram in [docs/ARCHITECTURE.md](ARCHITECTURE.md) and the layer-keyed
   anchors in [docs/CODEMAP.md](CODEMAP.md).
6. If visualization changed, regenerate example plots:

   ```bash
   uv run python examples/generate_example_plots.py
   ```

7. Push and open a PR against `main`. The CI workflow has four jobs — all must pass:
   - `lint` — `ruff check` and `ruff format --check` (Python 3.11, pinned ruff 0.15.2).
   - `docs-lint` — markdownlint, mermaid (maid + mmdc), yamllint, lychee link check.
   - `test` — pytest matrix across Python 3.10, 3.11, and 3.12; runs
     `uv run pytest -v -m "not integration"`.
   - `build` — `uv build` produces wheel and sdist artifacts.
8. There is no `.github/PULL_REQUEST_TEMPLATE.md` at time of writing — write a concise
   description covering *what* changed and *why*, and reference any related GitHub
   issues. Do not include AI-assistant attribution in commit messages or PR bodies.
9. Releases are cut from `main` by bumping `version` in `pyproject.toml`, running `uv lock`,
   changing the `## [Unreleased]` CHANGELOG heading to `## [X.Y.Z] - YYYY-MM-DD`, committing
   `pyproject.toml`, `uv.lock`, and `CHANGELOG.md` together, and creating a GitHub release
   with tag `vX.Y.Z`. `.github/workflows/publish.yml` then publishes to PyPI via Trusted
   Publishing. BiocondaBot opens a follow-up PR against bioconda-recipes automatically once
   the PyPI release is detected.

   Two things that catch people out:
   - `uv.lock` does not regenerate itself on a version bump. Skipping `uv lock` leaves a
     stale lock file committed alongside the new version.
   - Publishing triggers on the **GitHub release**, not on a tag push. Pushing a bare
     `vX.Y.Z` tag runs no workflow and publishes nothing.

# Development

Guidance for contributors working on pyLocusZoom locally. For first-run install and
quickstart, see [GETTING-STARTED.md](GETTING-STARTED.md) or the
[README](../README.md). For system architecture, see [ARCHITECTURE.md](ARCHITECTURE.md).
For environment variables and cache locations, see [CONFIGURATION.md](CONFIGURATION.md).

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
   uv run python -m pytest tests/
   ```

## Build Commands

pyLocusZoom has no `npm`-style `scripts` block; commands are run directly via `uv`.
The common development commands are:

| Command | Description |
|---------|-------------|
| `uv sync --all-extras` | Install all dependencies (runtime + `dev` + `spark` + `all` extras) into `.venv/`. |
| `uv run python -m pytest tests/` | Run the full test suite. Parallelism, the timeout, coverage and marker selection come from `addopts`. |
| `uv run python -m pytest tests/ --no-cov` | Fast iteration: skip coverage reporting. |
| `uv run python -m pytest tests/test_plotter_regional.py` | Run a single test file. |
| `uv run python -m pytest tests/ -m integration` | Run only the integration tests, which `addopts` deselects by default. |
| `uv tool run ruff check src/ tests/` | Run ruff lint checks (no fixes). |
| `uv tool run ruff format src/ tests/` | Apply ruff formatting. |
| `uv tool run ruff format --check src/ tests/` | Verify formatting without changes (matches CI check). |
| `uv run pre-commit run --all-files` | Run the full pre-commit suite against every file in the repo. |
| `uv build` | Build the wheel and sdist via hatchling into `dist/`. |
| `uv run python examples/generate_example_plots.py` | Regenerate example plots shown in the README. |
| `scripts/example_diff.sh [--keep]` | Generate outside the checkout and compare exports with HEAD. Exit 1 means differences; exit 2 means failure. `--keep` accepts generated changes only when affected exports have no manual edits. |
| `uv lock` | Refresh `uv.lock` after changing dependencies in `pyproject.toml`. |

See [CONTRIBUTING.md](../CONTRIBUTING.md) for the full pre-commit and pre-PR checklists.

## Example Exports

`examples/generate_example_plots.py` writes PNGs to `examples/matplotlib/` and HTML
to `examples/plotly/` and `examples/bokeh/`, relative to the working directory.
Both HTML backends load their JavaScript from a CDN (`include_plotlyjs="cdn"`,
bokeh's `CDN` resources), so no export embeds a runtime. Expect `examples/plotly/`
to total about 4.5 MB and `examples/bokeh/` about 5 MB; an export near 5 MB on its
own has the plotly.js bundle embedded. `scripts/example_diff.py` normalises
per-run div ids and treats a changed CDN version as a real difference, because it
changes the runtime the export loads.

The generator treats any `UserWarning` as an error. pyLocusZoom reports a skipped
layer (recombination overlay, gene track, LD colouring) as a `UserWarning` and
draws the figure without it, so a warning means a degraded export. The generator
fetches the canine recombination maps into the managed cache before plotting, so
the first run needs network; after that it runs offline. To fetch the maps ahead
of time:

```bash
uv run python -c "from pylocuszoom import download_canine_recombination_maps as d; print(d())"
```

`scripts/example_diff.sh` exits 2 when the generator fails, with or without
`--keep`, so a degraded export is never compared or accepted.

## Example Notebook

`examples/getting_started.ipynb` runs offline and seeds every cell that draws random
numbers, so a top-to-bottom run gives the same figures each time. The `dev`
dependency group provides `nbclient` and `ipykernel`. To execute the notebook
headless and fail on the first cell error:

```bash
uv run python scripts/execute_notebook.py examples/getting_started.ipynb
```

To execute it and save the outputs into the notebook, pass the same path as `--output`:

```bash
uv run python scripts/execute_notebook.py examples/getting_started.ipynb \
  --output examples/getting_started.ipynb
```

Commit saved outputs only when they come from a full top-to-bottom run like this one.

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

For a concrete example, see `prepare_manhattan_frames` in
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
| `pytest-cov` | Runs `uv run python -m pytest` on every Python change |

## Branch Conventions

- The default branch is **`main`**. CI runs on `push` to `main` and on all pull requests
  targeting `main` (`.github/workflows/ci.yml`).
- Feature branches are named `feature/your-feature` (`git checkout -b feature/your-feature`).
  Use `fix/…` for bug fixes and `chore/…` for maintenance tasks. No stricter convention
  is enforced. [CONTRIBUTING.md](../CONTRIBUTING.md#pull-request-guidelines) covers what a
  pull request needs.
- **Always commit and push before ending a session.** Do not leave `uv.lock` or other
  regenerated files uncommitted after a dependency bump.

## PR Process

1. Work on a feature branch off `main` (see above).
2. Write tests **before** the implementation where practical
   (see [docs/TESTING.md](TESTING.md) on test-driven development).
3. Before opening the PR, run the pre-commit checklist locally:

   ```bash
   uv run python -m pytest tests/
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

7. Push and open a PR against `main`. The CI workflow has five jobs — all must pass:
   - `lint` — `ruff check` and `ruff format --check` (Python 3.11, pinned ruff 0.15.2).
   - `docs-lint` — markdownlint, mermaid (maid + mmdc), yamllint, lychee link check.
   - `test` — pytest matrix across Python 3.10, 3.11, and 3.12; runs
     `uv run pytest`, which takes every flag from `addopts`.
   - `examples` — downloads the canine recombination maps, runs the example generator
     into a temporary directory (any `UserWarning` fails it) and executes the notebook.
   - `build` — `uv build` produces wheel and sdist artifacts.
8. There is no `.github/PULL_REQUEST_TEMPLATE.md` at time of writing — write a concise
   description covering *what* changed and *why*, and reference any related GitHub
   issues. Do not include AI-assistant attribution in commit messages or PR bodies.
9. Releases are cut from `main` by bumping `version` in `pyproject.toml`, running `uv lock`,
   changing the `## [Unreleased]` CHANGELOG heading to `## [X.Y.Z] - YYYY-MM-DD`, setting
   `version` and `date-released` in `CITATION.cff` to match, committing
   `pyproject.toml`, `uv.lock`, `CHANGELOG.md`, and `CITATION.cff` together, and creating a GitHub release
   with tag `vX.Y.Z`. `.github/workflows/publish.yml` then publishes to PyPI via Trusted
   Publishing. BiocondaBot opens a follow-up PR against bioconda-recipes automatically once
   the PyPI release is detected. Zenodo archives the release and mints a version DOI from
   `.zenodo.json`; the concept DOI in the README badge covers every version.

   Two things that catch people out:
   - `uv.lock` does not regenerate itself on a version bump. Skipping `uv lock` leaves a
     stale lock file committed alongside the new version.
   - Publishing triggers on the **GitHub release**, not on a tag push. Pushing a bare
     `vX.Y.Z` tag runs no workflow and publishes nothing.

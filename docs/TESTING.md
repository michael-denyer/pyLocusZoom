# Testing

This guide describes the test framework, how to run tests, and how to write new tests for pyLocusZoom.

## Test Framework and Setup

pyLocusZoom uses **pytest** (`>=7.0.0`) along with several plugins configured in `pyproject.toml` under `[project.optional-dependencies].dev`:

| Plugin | Version | Purpose |
|--------|---------|---------|
| `pytest` | `>=7.0.0` | Core test runner |
| `pytest-cov` | `>=4.0.0` | Coverage reporting (source: `pylocuszoom`, branch coverage enabled) |
| `pytest-randomly` | `>=3.0.0` | Randomizes test order to surface hidden ordering dependencies |
| `pytest-xdist` | `>=3.0.0` | Parallel test execution (`-n 3`) |
| `pytest-timeout` | `>=2.0.0` | Per-test timeout of 30s (catches hung tests and accidental network calls) |
| `hypothesis` | `>=6.0.0` | Property-based testing |

Shared fixtures (`sample_gwas_df`, `sample_genes_df`, `sample_exons_df`, `sample_recomb_df`) and Hypothesis profiles (`ci`, `dev`, `debug`) are defined in `tests/conftest.py`. The active Hypothesis profile is controlled by the `HYPOTHESIS_PROFILE` environment variable and defaults to `dev` (20 examples); CI typically sets `HYPOTHESIS_PROFILE=ci` (100 examples).

### Installing Test Dependencies

```bash
uv sync --extra dev --extra all
```

The `all` extra pulls in `pyspark` so PySpark-dependent tests can run.

## Running Tests

The `[tool.pytest.ini_options]` section of `pyproject.toml` sets default options:

```toml
addopts = "-n 3 --timeout=30 --cov=pylocuszoom --cov-report=term-missing -v -m 'not integration'"
```

That means `uv run pytest` already runs with three parallel workers, a 30-second per-test timeout, coverage reporting, verbose output, and integration tests deselected.

| Command | What It Runs |
|---------|--------------|
| `uv run python -m pytest tests/` | Full default suite (parallel, with coverage, integration tests skipped) |
| `uv run python -m pytest tests/ -n 3 --no-cov` | Fast iteration without coverage overhead |
| `uv run python -m pytest tests/test_plotter.py` | Single test file |
| `uv run python -m pytest tests/test_plotter.py::TestPlotEdgeCases -v` | Single test class |
| `uv run python -m pytest tests/test_plotter.py::TestPlotEdgeCases::test_name` | Single test function |
| `uv run python -m pytest -m integration` | Only integration tests (e.g. `tests/test_ensembl_integration.py`, which hits the live Ensembl API) |
| `uv run python -m pytest -p no:randomly` | Disable randomization to reproduce a specific order |
| `uv run python -m pytest --randomly-seed=<n>` | Reproduce a failing randomized run by replaying its seed |

### Markers

Currently a single custom marker is registered in `pyproject.toml`:

- `integration` — tests that require external services (e.g. the Ensembl REST API). Deselected by default; opt in with `-m integration`.

## Writing New Tests

### File Layout and Naming

Tests live under `tests/`. Files follow the `test_*.py` naming convention and map 1:1 onto source modules where possible. Examples:

- `tests/test_plotter.py` — regional `LocusZoomPlotter` tests
- `tests/test_manhattan_plotter.py`, `tests/test_qq.py`, `tests/test_manhattan.py` — Manhattan/QQ coverage
- `tests/test_stats_plotter.py`, `tests/test_phewas.py`, `tests/test_forest.py` — statistical plots
- `tests/test_ld.py` — PLINK wrapper (PLINK calls are mocked; no real PLINK binary required)
- `tests/test_notebook_backends.py` — Plotly/Bokeh backend parity checks
- `tests/test_recombination.py` — recombination map loading and CanFam4 liftover
- `tests/test_gene_track.py`, `tests/test_labels.py`, `tests/test_colors.py`, etc.

Test classes use `TestThing` and test functions use `test_behavior_when_condition`.

### Fixtures

Prefer the shared fixtures in `tests/conftest.py` (`sample_gwas_df`, `sample_genes_df`, `sample_exons_df`, `sample_recomb_df`) over constructing DataFrames inline. They use a seeded `numpy.random.default_rng(42)` so output is deterministic across randomized runs.

Hypothesis strategies shared across tests live in `tests/strategies.py`.

### Guidelines

- **Assert on observable outputs, not mock call counts.** Check returned figures, DataFrame columns/shapes, written files, and raised exceptions. Reserve `assert_called_once_with` for true system boundaries (PLINK subprocess, HTTP, filesystem dispatch).
- **Mock PLINK calls** — tests must not require a real PLINK installation. See `tests/test_ld.py` for the established pattern.
- **Cover edge cases**: empty DataFrames, missing required columns, mismatched list lengths, single-SNP regions, and cross-chromosome filtering.
- **Respect the 30s timeout.** If a test is legitimately slow, override with `@pytest.mark.timeout(60)` rather than raising the global default.
- **Randomization-safe**: tests must not depend on execution order. If a test only passes under a specific seed, that is a bug in the test.

## Coverage Requirements

Coverage is measured with `pytest-cov` over the `pylocuszoom` package with branch coverage enabled (`[tool.coverage.run]` in `pyproject.toml`). The terminal report lists missing lines via `--cov-report=term-missing`.

The following patterns are excluded from coverage:

- `pragma: no cover`
- `if TYPE_CHECKING:`
- `raise NotImplementedError`

**No minimum coverage threshold is configured.** There is no `fail_under` in `[tool.coverage.report]` and CI does not enforce a numeric floor; coverage is reported for visibility only.

## CI Integration

Tests run in the `test` job of `.github/workflows/ci.yml`, triggered on pushes and pull requests targeting `main`. The job uses a matrix across Python 3.10, 3.11, and 3.12 on `ubuntu-latest`.

Steps:

1. Check out the repo.
2. Install `uv` via `astral-sh/setup-uv`.
3. `uv python install ${{ matrix.python-version }}`.
4. `uv sync --extra dev --extra all` to install dev and PySpark dependencies.
5. `uv run pytest -v -m "not integration"` to run the suite (integration tests are skipped in CI because they hit the live Ensembl API).

Separate jobs in the same workflow handle linting (`ruff check`, `ruff format --check` pinned to `ruff@0.15.2`), documentation linting (markdownlint, mermaid maid + renderer parity, yamllint, lychee link check), and package building (`uv build`). A test failure, lint failure, or doc-lint failure will block the PR.

Because `pytest-xdist` and `pytest-randomly` are active, every CI run reports the worker count and the random seed in the header — use `pytest --randomly-seed=<seed>` locally to reproduce a failure.

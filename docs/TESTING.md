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

`tests/test_rendering_contract.py` drives the plotters through a `RecordingBackend` and asserts on the call sequence it records, because no backend can serialise a figure the same way twice. Every other test asserts on observable output. A change to a backend or a panel must also pass `scripts/example_diff.sh`. It regenerates `examples/`, normalises plotly UUIDs and bokeh element ids, and prints `NO REAL DIFFS` or one `REAL DIFF:` line per export whose content changed. Matplotlib PNGs are deterministic, so any PNG diff is a real rendering change.

Shared fixtures and Hypothesis profiles (`ci`, `dev`, `debug`) are defined in `tests/conftest.py`. The active Hypothesis profile is controlled by the `HYPOTHESIS_PROFILE` environment variable and defaults to `dev` (20 examples). CI sets `HYPOTHESIS_PROFILE=ci` (100 examples), which costs under a second on the full suite.

### Installing Test Dependencies

```bash
uv sync --extra dev --extra all
```

The `all` extra pulls in `pyspark` so PySpark-dependent tests can run.

## Running Tests

Defaults come from `[tool.pytest.ini_options].addopts` in
[`pyproject.toml`](../pyproject.toml), which is the only place they are
written down. `uv run pytest` already runs with parallel workers, a per-test
timeout, coverage reporting, verbose output, and integration tests deselected,
so none of those flags belong on a command line or in a CI step.

| Command | What It Runs |
|---------|--------------|
| `uv run python -m pytest tests/` | Full default suite (parallel, with coverage, integration tests skipped) |
| `uv run python -m pytest tests/ --no-cov` | Fast iteration without coverage overhead |
| `uv run python -m pytest tests/test_plotter.py` | Single test file |
| `uv run python -m pytest tests/test_plotter.py::TestPlotEdgeCases` | Single test class |
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

- `tests/test_plotter.py` — the regional `LocusZoomPlotter` API
- `tests/test_plotter_ld.py` — LD calculation and the LD heatmap panel
- `tests/test_plotter_recombination.py` — the recombination overlay and its download errors
- `tests/test_plotter_backends.py` — regional plots rendered through each backend
- `tests/test_data_intake.py` — the shared p-value intake policy
- `tests/test_manhattan_plotter.py`, `tests/test_qq.py`, `tests/test_manhattan.py` — Manhattan/QQ coverage
- `tests/test_stats_plotter.py`, `tests/test_phewas.py`, `tests/test_forest.py` — statistical plots
- `tests/test_ld.py` — PLINK wrapper (driven through the `fake_plink` fixture; no real PLINK binary required)
- `tests/test_backends.py` — the shared `PlotBackend` surface and the matplotlib backend
- `tests/test_plotly_backend.py`, `tests/test_bokeh_backend.py` — the interactive backends
- `tests/test_notebook_backends.py` — Plotly/Bokeh notebook compatibility
- `tests/test_recombination.py` — recombination map loading and CanFam4 liftover
- `tests/test_ucsc.py` — the UCSC gene client and the build-to-source routing (HTTP mocked)
- `tests/test_fixture_hygiene.py` — fails if one fixture name gains a second schema
- `tests/test_docs_contract.py` — fails if a documented pytest command repeats an `addopts` flag, or if `__all__` outgrows the CODEMAP table
- `tests/test_coerce.py`, `tests/test_plotly_layout.py` — the pure helpers the interactive backends share
- `tests/test_ld_plotting.py` — the policy deciding whether a plot reaches for PLINK
- `tests/test_loaders.py` — loader dispatch, format detection and file-path validation
- `tests/test_loaders_gwas.py`, `tests/test_loaders_eqtl.py`, `tests/test_loaders_finemapping.py`, `tests/test_loaders_annotation.py` — one file per loader family
- `tests/test_validation_contract.py` — the load-time contract as one accepted table and one rejected table; a new column rule is a new row
- `tests/test_validation.py` — the `ColumnSpec` rule engine the contract runs on
- `tests/test_gene_track.py`, `tests/test_labels.py`, `tests/test_colors.py`, etc.

Do not name a file after a batch of bugs. A maintainer editing `bokeh_backend.py`
must be able to find its tests from the module name alone.

Test classes use `TestThing` and test functions use `test_behavior_when_condition`.

### Fixtures

Prefer the shared fixtures in `tests/conftest.py` over constructing DataFrames inline. They use a seeded `numpy.random.default_rng(42)` so output is deterministic across randomized runs.

One fixture name means exactly one schema. `regional_gwas_df`, `small_regional_gwas_df` and `tiny_regional_gwas_df` are the `rs`/`pos`/`p_value` region shapes; `manhattan_gwas_df` and `manhattan_rs_gwas_df` are the `chr`/`pos`/`p_value` shapes; `labelled_gwas_df` carries a precomputed `neglog10p`. `test_fixture_hygiene.py` fails if any name gains a second schema, so a new shape needs a new name rather than a local shadow. It reads the column names of every `pd.DataFrame({...})` literal in the fixture body, so a fixture that assembles several frames and returns a list is covered too.

An autouse fixture closes every pyplot figure after each test, so a test that
builds a figure does not need to close it.

`warning_records` collects `pylocuszoom` warnings. loguru does not feed pytest's `caplog`, so a test that takes `caplog` and asserts on it will pass no matter what the code logs.

Hypothesis strategies shared across tests live in `tests/strategies.py`.

### Guidelines

- **Assert on observable outputs, not mock call counts.** Check returned figures, DataFrame columns/shapes, written files, and raised exceptions. Reserve `assert_called_once_with` for true system boundaries (PLINK subprocess, HTTP, filesystem dispatch).
- **Drive PLINK through `fake_plink`** — tests must not require a real PLINK installation. The `fake_plink` fixture in `conftest.py` patches `subprocess.run` and writes a real `.ld` file at the path the command asked for, so command construction, output parsing and the R2 merge all stay inside the test. Assert on the frame `calculate_ld` returns, not on what the mock received: a command flag is already pinned by `TestBuildLdCommand` and `TestBuildPairwiseLdCommand`, which call the pure builders and assert on the list they return.
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
5. `uv run pytest` to run the suite. Every flag comes from `addopts`, including the marker expression that deselects the integration tests, which hit the live Ensembl API.

Separate jobs in the same workflow handle linting (`ruff check`, `ruff format --check` pinned to `ruff@0.15.2`), documentation linting (markdownlint, mermaid maid + renderer parity, yamllint, lychee link check), and package building (`uv build`). A test failure, lint failure, or doc-lint failure will block the PR.

Because `pytest-xdist` and `pytest-randomly` are active, every CI run reports the worker count and the random seed in the header — use `pytest --randomly-seed=<seed>` locally to reproduce a failure.

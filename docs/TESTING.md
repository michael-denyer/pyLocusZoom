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
| `pytest-timeout` | `>=2.0.0` | Per-test timeout of 30s (catches hung tests) |
| `hypothesis` | `>=6.0.0` | Property-based testing |

`tests/test_rendering_contract.py` drives the plotters through a `RecordingBackend` and asserts on the call sequence it records, because no backend can serialise a figure the same way twice. Every other test asserts on observable output. A change to a backend or a panel must also pass `scripts/example_diff.sh`. It regenerates `examples/`, normalises plotly UUIDs and bokeh element ids, and prints `NO REAL DIFFS` or one `REAL DIFF:` line per export whose content changed. Matplotlib PNGs are deterministic, so any PNG diff is a real rendering change.

Shared fixtures and Hypothesis profiles (`ci`, `dev`, `debug`) are defined in `tests/conftest.py`. The active Hypothesis profile is controlled by the `HYPOTHESIS_PROFILE` environment variable and defaults to `dev` (20 examples). CI sets `HYPOTHESIS_PROFILE=ci` (100 examples), which costs under a second on the full suite.

### Installing Test Dependencies

```bash
uv sync --extra dev --extra all
```

The `dev` extra is all the suite needs: the PySpark tests use a stand-in frame
that answers `toPandas()`, so `uv sync --extra dev` runs every test. CI also
installs the `all` extra, which adds `pyspark`.

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
| `uv run python -m pytest tests/test_plotter_regional.py` | Single test file |
| `uv run python -m pytest tests/test_plotter_regional.py::TestLeadPosBoundary` | Single test class |
| `uv run python -m pytest tests/test_plotter_regional.py::TestLeadPosBoundary::test_name` | Single test function |
| `uv run python -m pytest -m integration` | Only integration tests (e.g. `tests/test_ensembl_integration.py`, which hits the live Ensembl API) |
| `uv run python -m pytest -p no:randomly` | Disable randomization to reproduce a specific order |
| `uv run python -m pytest --randomly-seed=<n>` | Reproduce a failing randomized run by replaying its seed |

### Markers

Currently a single custom marker is registered in `pyproject.toml`:

- `integration` — tests that require external services (e.g. the Ensembl REST API). Deselected by default; opt in with `-m integration`. Every other test runs with outbound connections and DNS lookups refused by the autouse `block_network` fixture in `tests/conftest.py`, so a test that would download fails instead of reading or filling the machine's cache.

## Writing New Tests

### File Layout and Naming

Tests live under `tests/`. Files follow the `test_*.py` naming convention and map 1:1 onto source modules where possible.

**One subject, one file.** Every public subject has exactly one owning file, and a test for that subject goes there, not into whichever file was open. Before adding a test, find the owner below; if a behaviour already has a test there, strengthen that test rather than writing a second one elsewhere. A subject may span two files only when one would pass 1,000 lines, and then the split is by sub-subject, named in this list. No test file exceeds 1,000 lines.

| Subject | Owning file |
|---|---|
| Regional `LocusZoomPlotter`: construction, plot options, region selection, lead resolution, stacked-list validation | `tests/test_plotter_regional.py` |
| Automatic gene fetching for the regional gene track | `tests/test_plotter_genes.py` |
| LD colouring and the LD heatmap panel in regional plots | `tests/test_plotter_ld.py` |
| The recombination overlay in regional plots | `tests/test_plotter_recombination.py` |
| Regional plots rendered through each backend | `tests/test_plotter_backends.py` |
| Each regional panel type building and drawing itself | `tests/test_regional_plan.py` |
| `render_figure` and the figure-level primitives a `FigurePlan` names | `tests/test_figure.py` |
| Every public plot method collecting its frames through `to_pandas` | `tests/test_spark_intake.py` |
| Threshold defaults, `UNSET` and `add_significance_line` | `tests/test_plotter_utils.py` |
| Shared p-value intake: NaN, out-of-range, empty and missing-column input | `tests/test_data_intake.py` |
| Config models, including the regional option surface | `tests/test_config.py` |
| `ManhattanPlotter`: construction, input boundary, threshold defaults | `tests/test_manhattan_plotter.py` |
| `ManhattanPlotter`: output of each `plot_*` method | `tests/test_manhattan_plotter_methods.py` |
| Manhattan and QQ pure data preparation | `tests/test_manhattan.py`, `tests/test_qq.py` |
| `GenomeWideStyle` across the genome-wide plotters | `tests/test_genomewide_style.py` |
| `MiamiPlotter` | `tests/test_miami_plotter.py` |
| `ColocPlotter` | `tests/test_coloc_plotter.py`; the `coloc_plot_spec` contract in `tests/test_coloc.py` |
| `LDHeatmapPlotter` | `tests/test_ld_heatmap_plotter.py` |
| `StatsPlotter`, PheWAS, forest | `tests/test_stats_plotter.py`, `tests/test_phewas.py`, `tests/test_forest.py` |
| Fine-mapping data and the fine-mapping panel | `tests/test_finemapping.py` |
| eQTL filtering and `calculate_colocalization_overlap` | `tests/test_eqtl.py` |
| Gene row layout, strand arrows and the gene panel | `tests/test_gene_track.py` |
| SNP label ranking and placement | `tests/test_labels.py` |
| The palette in `colors.py`, and that no other module holds a hex literal | `tests/test_colors.py` |
| The `PlotBackend` surface shared by all backends, and the matplotlib backend | `tests/test_backends.py` |
| Plotly-only and bokeh-only library behaviour | `tests/test_plotly_backend.py`, `tests/test_bokeh_backend.py` |
| Legend, recombination-overlay and heatmap composition above the seam | `tests/test_composition.py` |
| Hover columns, roles and each backend's tooltip spec | `tests/test_hover.py` |
| Coercions out of matplotlib's vocabulary; plotly subplot geometry | `tests/test_coerce.py`; `tests/test_plotly_layout.py` |
| Notebook and HTML export of the interactive backends | `tests/test_notebook_backends.py` |
| Rendering call sequence through `RecordingBackend` | `tests/test_rendering_contract.py` |
| PLINK command construction, output parsing, process execution | `tests/test_ld.py`, `tests/test_ld_parsing.py`, `tests/test_ld_process.py` |
| Whether a plot reaches for PLINK | `tests/test_ld_enrichment.py` |
| Recombination map loading, region lookup and liftover | `tests/test_recombination.py` |
| Fetching, unpacking and publishing the managed map set | `tests/test_recombination_maps.py` |
| Coordinate liftover | `tests/test_liftover.py` |
| Species records and `resolve_species` | `tests/test_species.py` |
| Ensembl gene source; the same against the live API (`integration`) | `tests/test_ensembl.py`; `tests/test_ensembl_integration.py` |
| UCSC gene source and the build-to-source router | `tests/test_ucsc.py` |
| The gene and exon cache | `tests/test_gene_cache.py` |
| The retrying HTTP transport | `tests/test_http.py` |
| Genome-build records and name folding | `tests/test_genome_build.py` |
| Exception hierarchy | `tests/test_exceptions.py` |
| Package exports | `tests/test_init.py` |
| The canonical column vocabulary | `tests/test_canonical_columns.py` |
| `utils`: `to_pandas`, chromosome normalising, region filtering, the cache base | `tests/test_utils.py` |
| Logging switches and sinks | `tests/test_logging.py` |
| Loaders | `tests/test_loaders.py` (dispatch, format detection, file paths), one `tests/test_loaders_<family>.py` per loader family |
| Load-time column contract; the rule engine it runs on | `tests/test_validation_contract.py`; `tests/test_validation.py` |
| `scripts/example_diff.sh` | `tests/test_example_diff_script.py` |
| Suite structure: fixture schemas; documentation contracts (documented commands, the CODEMAP and USER_GUIDE public-API tables) | `tests/test_fixture_hygiene.py`; `tests/test_docs_contract.py` |

`tests/figure_probes.py` is the one probe object per backend (`PROBES`). It translates panel count, tick labels, y-axis range, legend corner and swatch edges, horizontal and vertical lines, rectangles, region highlights, colour-bar titles, scatter marker positions, point alpha, font sizes, marker symbols and hover (values, and each tooltip's fields and number formats) into one vocabulary, and it is the only place in the suite that knows matplotlib's, plotly's or bokeh's figure internals. `marker_symbols`, `has_hover`, `hover_values`, `hover_fields`, `standalone_html` and `json_payload` exist only for the interactive backends. A matplotlib-only test may read the matplotlib `Figure` directly.

### Private seams tests may touch

Tests assert on public behaviour. These private names are the deliberate exceptions, each because the behaviour has no public handle or because the ADRs make the seam part of the design. A new private import or patch target outside this list needs a reason in review, and ideally a public handle instead.

- **Plan layer** ([ADR 0001](adr/0001-deepen-rendering-seam.md), [ADR 0007](adr/0007-one-figure-plan.md)): `pylocuszoom._figure` (`FigurePlan`, `render_figure`, `RegionHighlight`), the panel classes in `pylocuszoom.panels.*`, built and drawn directly in `test_regional_plan.py`, `test_finemapping.py` and `test_gene_track.py` (including `MiamiRequest`, `AssociationInput` and `AssociationPanel`), and the `RecordingBackend` in `test_rendering_contract.py`.
- **Deep internal modules with their own unit tests**, private helpers included (such as `_gene_cache._entry_file` and `_http._with_retries`, and `_http._stream_to` as a patch target): `_gene_cache`, `_gene_source`, `_http`, `_data`, `_liftover`, `_ld_enrichment`, `_plotter_utils`, `panels._shared` (`add_significance_line`, in `tests/test_plotter_utils.py`), `backends._coerce`, `plotly_layout`.
- **I/O boundaries patched where they are looked up**: `pylocuszoom._http.requests.get` and `_http.time.sleep` (HTTP), `subprocess.run` through `fake_plink` (PLINK), `pylocuszoom.ld.find_plink`, `pylocuszoom._ld_enrichment.calculate_ld`, the `pylocuszoom.recombination` download and directory functions, and `pylocuszoom._liftover.download_file` (chain downloads). A test that reads or fills the managed cache takes the `cache_home` fixture, which points `XDG_CACHE_HOME` at the test's directory.
- **Pure private helpers pinned for their edge cases**: `recombination._stage_archive` and `_publish_map_generation` (archive safety), `ld._resolve_plink` and `_add_species_flags`, `loaders.gwas._detect_format`, `colors._find_eqtl_bin`, `panels.coloc._get_effect_agreement_color`, `bokeh_backend._create_color_palette`, `utils._platform_cache_base`.
- **Registry and logger state**: `backends._BACKENDS` (registering a test backend), and loguru's `logger.enable`/`logger.disable("pylocuszoom")` and `logger._core.handlers` (switching the package's records on for `conftest.warning_records`, and checking that import leaves the handlers alone).

`LocusZoomPlotter` private methods and caches, and the coloc merge and lead functions, are not seams: test them through the rendered figure. A panel class is a seam only for building and drawing that panel on its own; what a plot shows the user is tested through the rendered figure.

Do not name a file after a batch of bugs. A maintainer editing `bokeh_backend.py`
must be able to find its tests from the module name alone.

Test classes use `TestThing` and test functions use `test_behavior_when_condition`.

A test's name is a claim about what it checks. `test_significance_lines` must assert a line was drawn; if the body only proves the call returned, the name lies and a reviewer reading the test list believes the behaviour is pinned when it is not. Where the behaviour really is "this input does not crash", say so in the name (`test_..._does_not_raise`) and drop the assertion rather than asserting the figure exists.

### Fixtures

Prefer the shared fixtures in `tests/conftest.py` over constructing DataFrames inline. They use a seeded `numpy.random.default_rng(42)` so output is deterministic across randomized runs.

One fixture name means exactly one schema. `regional_gwas_df`, `small_regional_gwas_df` and `tiny_regional_gwas_df` are the `rs`/`pos`/`p_value` region shapes; `manhattan_gwas_df` and `manhattan_rs_gwas_df` are the `chr`/`pos`/`p_value` shapes; `labelled_gwas_df` carries a precomputed `neglog10p`. `test_fixture_hygiene.py` fails if any name gains a second schema, so a new shape needs a new name rather than a local shadow. It reads the column names of every `pd.DataFrame({...})` literal in the fixture body, so a fixture that assembles several frames and returns a list is covered too.

An autouse fixture closes every pyplot figure after each test, so a test that
builds a figure does not need to close it.

`warning_records` collects `pylocuszoom` WARNING log records, switching the package's logging on for the test and off again. loguru does not feed pytest's `caplog`, so a test that takes `caplog` and asserts on it will pass no matter what the code logs. A layer the plot skips is a `UserWarning`, not a log record: assert it with `pytest.warns`.

Hypothesis strategies shared across tests live in `tests/strategies.py`.

### Guidelines

- **Assert on observable outputs, not mock call counts.** Check returned figures, DataFrame columns/shapes, written files, and raised exceptions. Reserve `assert_called_once_with` for true system boundaries (PLINK subprocess, HTTP, filesystem dispatch).
- **Drive PLINK through `fake_plink`** — tests must not require a real PLINK installation. The `fake_plink` fixture in `conftest.py` patches `subprocess.run` and writes a real `.ld` file at the path the command asked for, so command construction, output parsing and the R2 assignment all stay inside the test. Assert on the frame `calculate_ld` returns, not on what the mock received: a command flag is already pinned by `TestBuildLdCommand` and `TestBuildPairwiseLdCommand`, which call the pure builders and assert on the list they return.
- **State a rendering behaviour once, not once per backend.** A fact about the figure (how many panels, where the threshold line sits, which marker, whether it hovers) belongs in one `@pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)` test reading a probe from `tests/figure_probes.py`; use `INTERACTIVE_BACKENDS` only for questions matplotlib cannot answer, such as HTML export. If the probe cannot answer the question yet, add a method to all three probes rather than walking the figure in the test. Hand-written per-backend twins drift: the bokeh eQTL marker test used to pass with the negative-effect glyph never drawn. A genuine library-specific regression still belongs in `test_plotly_backend.py` or `test_bokeh_backend.py`.
- **Cover edge cases**: empty DataFrames, missing required columns, mismatched list lengths, single-SNP regions, and cross-chromosome filtering.
- **Respect the 30s timeout.** If a test is legitimately slow, override with `@pytest.mark.timeout(60)` rather than raising the global default.
- **Randomization-safe**: tests must not depend on execution order. If a test only passes under a specific seed, that is a bug in the test.

## Model Checking

`specs/tla/RecombPublish.tla` models concurrent writers and a reader of the recombination map cache, one filesystem call per step, and checks the claims in the `_publish_map_generation` docstring. It is not run in CI. After changing `_publish_map_generation` or its callers, update the spec to match and run both configurations with TLC ([`tla2tools.jar`](https://github.com/tlaplus/tlaplus/releases) and Java 11 or later):

```bash
cd specs/tla
java -cp tla2tools.jar tlc2.TLC -config RecombPublish.cfg RecombPublish.tla
java -cp tla2tools.jar tlc2.TLC -config RecombPublishNoGap.cfg RecombPublish.tla
```

Each run ends with `Model checking completed. No error has been found.` A violation prints the interleaving that breaks the claim; replay it against the real function with a test that patches the check the model split from its action, as `test_a_writer_that_loses_the_legacy_symlink_race_still_succeeds` does.

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

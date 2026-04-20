# Configuration

pyLocusZoom is a Python library, so most "configuration" is done through
keyword arguments passed to the plotting API rather than through environment
variables or config files. This document covers the three places where
behavior can be configured:

1. **Environment variables** that control on-disk cache locations.
2. **Programmatic configuration** via the Pydantic models in
   [`src/pylocuszoom/config.py`](../src/pylocuszoom/config.py).
3. **Project-level configuration** declared in
   [`pyproject.toml`](../pyproject.toml) (for developers working on the
   package itself).

There is no `.env` file and no runtime configuration file. The library does
not read any `PYLOCUSZOOM_*` environment variables.

## Environment Variables

pyLocusZoom honours a small number of standard OS-level environment
variables that control where cached reference data (recombination maps,
Ensembl gene annotations) is stored. These are the only environment
variables the library reads.

| Variable          | Required | Default                                                                                  | Description                                                                                                                                                 |
| ----------------- | -------- | ---------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `XDG_CACHE_HOME`  | Optional | `~/.cache` (Linux) / `~/.cache` (macOS)                                                  | Base directory for the cache on Linux. When set, caches are placed under `$XDG_CACHE_HOME/pylocuszoom/`. Not consulted on macOS by `ensembl.get_ensembl_cache_dir()`. |
| `LOCALAPPDATA`    | Optional | `%USERPROFILE%` (recombination) / `%USERPROFILE%\AppData\Local` (ensembl)                | Base directory for the cache on Windows. Caches are placed under `%LOCALAPPDATA%\pylocuszoom\`.                                                             |

Resolved cache paths by platform:

| Platform   | Recombination maps                                                                        | Ensembl gene cache                                       |
| ---------- | ----------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| macOS      | `~/.cache/pylocuszoom/recombination_maps` (or `$XDG_CACHE_HOME/pylocuszoom/recombination_maps`) | `~/.cache/pylocuszoom/ensembl`                           |
| Linux      | `$XDG_CACHE_HOME/pylocuszoom/recombination_maps` or `~/.cache/pylocuszoom/recombination_maps` | `$XDG_CACHE_HOME/pylocuszoom/ensembl` or `~/.cache/pylocuszoom/ensembl` |
| Windows    | `%LOCALAPPDATA%\pylocuszoom\recombination_maps`                                           | `%LOCALAPPDATA%\pylocuszoom\ensembl`                     |
| Databricks | `/dbfs/FileStore/reference_data/recombination_maps` (auto-detected when `/dbfs` exists)   | Falls through to the Linux/XDG path                      |

Implementation:

- [`recombination.get_default_data_dir()`](../src/pylocuszoom/recombination.py)
- [`ensembl.get_ensembl_cache_dir()`](../src/pylocuszoom/ensembl.py)

You can also override the cache location explicitly by passing `output_dir`
to `download_canine_recombination_maps()` / `ensure_recomb_maps()` — this
bypasses the environment variables entirely.

## Programmatic Configuration (Pydantic Models)

The user-facing API uses plain keyword arguments (`plot()`,
`plot_stacked()`). Internally these kwargs are validated by frozen
Pydantic models defined in
[`src/pylocuszoom/config.py`](../src/pylocuszoom/config.py). You normally
do not construct these directly, but they define the canonical set of
options and their defaults.

### `RegionConfig` — genomic region (required)

| Field   | Type          | Default  | Validation                                     |
| ------- | ------------- | -------- | ---------------------------------------------- |
| `chrom` | `int \| str`  | required | Integer `>= 1`, or non-empty string            |
| `start` | `int`         | required | `>= 0`                                         |
| `end`   | `int`         | required | `> 0` and strictly greater than `start`        |

### `ColumnConfig` — GWAS DataFrame column names

| Field     | Type  | Default   | Description           |
| --------- | ----- | --------- | --------------------- |
| `pos_col` | `str` | `"ps"`    | Position column name  |
| `p_col`   | `str` | `"p_wald"`| P-value column name   |
| `rs_col`  | `str` | `"rs"`    | SNP identifier column |

### `DisplayConfig` — visual options

| Field                | Type                     | Default       | Description                           |
| -------------------- | ------------------------ | ------------- | ------------------------------------- |
| `snp_labels`         | `bool`                   | `True`        | Draw SNP labels                       |
| `label_top_n`        | `int` (`>= 0`)           | `5`           | Number of top SNPs to label           |
| `show_recombination` | `bool`                   | `True`        | Overlay recombination rate track      |
| `figsize`            | `tuple[float, float]`    | `(12.0, 8.0)` | Figure size in inches (width, height) |

### `LDConfig` — linkage disequilibrium

| Field               | Type              | Default | Description                                |
| ------------------- | ----------------- | ------- | ------------------------------------------ |
| `lead_pos`          | `int \| None`     | `None`  | Position of lead SNP (`>= 1` when set)     |
| `ld_reference_file` | `str \| None`     | `None`  | Path to PLINK binary fileset               |
| `ld_col`            | `str \| None`     | `None`  | Column with pre-computed R² values         |

Cross-field rules:

- `ld_col` and `ld_reference_file` are mutually exclusive.
- If `ld_reference_file` is set, `lead_pos` is required (enforced on
  `PlotConfig`).

### Composite configs

- `PlotConfig` composes `RegionConfig`, `ColumnConfig`, `DisplayConfig`, and
  `LDConfig`. Use `PlotConfig.from_kwargs(...)` to construct one from the
  same flat kwargs that `plot()` accepts.
- `StackedPlotConfig` extends the pattern with list-valued
  `lead_positions`, `panel_labels`, and `ld_reference_files` fields for
  multi-panel plots.
- `ColocConfig` covers colocalisation-specific options.

All config models are `frozen=True` — construct a new instance rather than
mutating an existing one.

## Required vs Optional Settings

Because configuration is passed at call time, "required" here means
"must be supplied when calling `plot()` / `plot_stacked()`".

| Setting                  | Required?                              | Notes                                                      |
| ------------------------ | -------------------------------------- | ---------------------------------------------------------- |
| `chrom`, `start`, `end`  | Required                               | Validation error if missing or if `start >= end`.          |
| `pos_col`, `p_col`, `rs_col` | Optional                           | Default to `"ps"`, `"p_wald"`, `"rs"`.                     |
| `lead_pos`               | Required *if* `ld_reference_file` set  | Otherwise optional.                                        |
| `ld_reference_file`      | Optional                               | Mutually exclusive with `ld_col`.                          |
| `ld_col`                 | Optional                               | Mutually exclusive with `ld_reference_file`.               |
| `snp_labels`, `label_top_n`, `show_recombination`, `figsize` | Optional | Sensible defaults (see table above).      |

Validation failures raise `pydantic.ValidationError` at call time.

## Defaults Summary

Defaults defined in source (see
[`config.py`](../src/pylocuszoom/config.py)):

```text
pos_col            = "ps"
p_col              = "p_wald"
rs_col             = "rs"
snp_labels         = True
label_top_n        = 5
show_recombination = True
figsize            = (12.0, 8.0)
lead_pos           = None
ld_reference_file  = None
ld_col             = None
```

## Per-Environment Overrides

pyLocusZoom does not distinguish "development" vs "production" environments
at runtime — it is a library, not a service. There are no
`.env.development` / `.env.production` files and no `NODE_ENV`-style switch.

If you need per-environment behaviour, do it at the caller level, e.g.:

- Set `XDG_CACHE_HOME` / `LOCALAPPDATA` per machine to control where
  reference data is cached.
- Pre-download reference data in CI with `ensure_recomb_maps()` pointing at
  a shared directory, then set `output_dir=` accordingly at runtime.
- On Databricks, the `/dbfs/FileStore/reference_data/recombination_maps`
  path is selected automatically.

## Project / Developer Configuration

The following settings live in [`pyproject.toml`](../pyproject.toml) and
only affect contributors working on pyLocusZoom itself (not library users):

| Setting                           | Value                                                                   |
| --------------------------------- | ----------------------------------------------------------------------- |
| `requires-python`                 | `>= 3.10`                                                               |
| `build-system.requires`           | `hatchling==1.29.0`                                                     |
| `tool.pytest.ini_options.addopts` | `-n 3 --timeout=30 --cov=pylocuszoom --cov-report=term-missing -v -m 'not integration'` |
| `tool.ruff.line-length`           | `88`                                                                    |
| `tool.ruff.target-version`        | `py310`                                                                 |
| `tool.ruff.lint.select`           | `["E", "F", "I", "W"]`                                                  |
| `tool.ruff.lint.ignore`           | `["E501"]`                                                              |

Optional dependency groups: `dev`, `spark`, `all` (see `pyproject.toml`).

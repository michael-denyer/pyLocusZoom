# Getting Started

This guide walks you through installing pyLocusZoom, producing a first plot, and troubleshooting common setup issues. pyLocusZoom is a publication-ready regional association plot library for GWAS data with LD coloring, gene tracks, and recombination overlays.

## Prerequisites

| Tool | Version | Why |
|------|---------|-----|
| Python | >=3.10 | Declared in `pyproject.toml` (`requires-python = ">=3.10"`); tested against 3.10, 3.11, 3.12 |
| `uv` (recommended) | latest | Project uses `uv.lock` for reproducible builds |
| `pip` (alternative) | any recent | Works for consumers who only want to install from PyPI |
| PLINK 1.9 | any | Optional — only needed if you compute LD from genotype files instead of passing a pre-computed `ld_col` |
| Node.js | >=20 | Only required for contributors — pre-commit hooks (`mermaid-lint`, `mermaid-render`) invoke `npx` |

All Python runtime dependencies (matplotlib, pandas, numpy, loguru, pyliftover, plotly, bokeh, kaleido, adjustText, pydantic, requests, tqdm, colorcet) are installed automatically.

## Installation

### As a user (from PyPI)

```bash
pip install pylocuszoom
```

Or with uv:

```bash
uv add pylocuszoom
```

Or with conda (Bioconda):

```bash
conda install -c bioconda pylocuszoom
```

### As a contributor (from source)

1. Clone the repository:

   ```bash
   git clone https://github.com/michael-denyer/pyLocusZoom.git
   cd pyLocusZoom
   ```

2. Install all dependencies (runtime + dev extras) with uv:

   ```bash
   uv sync --all-extras
   ```

   This installs the `dev` extras (`pytest`, `pytest-cov`, `pytest-randomly`, `pytest-xdist`, `pytest-timeout`, `hypothesis`, `ruff`) defined in `pyproject.toml`.

3. (Optional, contributors only) Install pre-commit hooks:

   ```bash
   uv run pre-commit install
   ```

## First Run

The fastest path from install to a working plot is a regional association plot against a small GWAS DataFrame.

```python
import pandas as pd
from pylocuszoom import LDConfig, LocusZoomPlotter

# Minimal GWAS DataFrame: chrom, pos, p_value (rs optional but recommended for labels)
gwas_df = pd.DataFrame({
    "chrom": [1] * 5,
    "pos":   [1_100_000, 1_250_000, 1_500_000, 1_750_000, 1_900_000],
    "p_value": [1e-3, 1e-5, 1e-8, 1e-6, 1e-4],
    "rs":    ["rs1", "rs2", "rs3", "rs4", "rs5"],
})

plotter = LocusZoomPlotter(species="canine", auto_genes=False)
fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=1_000_000,
    end=2_000_000,
    ld=LDConfig(lead_pos=1_500_000),
)
fig.savefig("regional_plot.png", dpi=150)
```

You can also explore the bundled Jupyter tutorial:

```bash
uv run jupyter notebook examples/getting_started.ipynb
```

Or regenerate all example plots (useful to confirm the install works end-to-end):

```bash
uv run python examples/generate_example_plots.py
```

Generated PNGs land in `examples/matplotlib/`, `examples/plotly/`, and `examples/bokeh/`.

## Common Setup Issues

- **`ModuleNotFoundError: No module named 'pylocuszoom'`** — You installed into a different environment than the one running Python. With uv, always prefix commands with `uv run` (e.g., `uv run python script.py`) so they execute inside the project's `.venv`.
- **Python version too old** — `pip install pylocuszoom` fails with a version-compatibility message on Python <3.10. Upgrade Python or use `uv python install 3.12`.
- **PLINK not found when computing LD** — Passing `ld_reference_file="genotypes"` requires PLINK 1.9 on your `PATH`, or an explicit `plink_path="/path/to/plink"`. If you already have LD values, pass them via the `ld_col` parameter instead and skip PLINK entirely.
- **First canine plot is slow / appears to hang** — The first regional plot for a canine locus downloads ~50MB of recombination maps into the user cache. Subsequent plots are instant. You can pre-download via `pylocuszoom.recombination.ensure_recomb_maps()`.
- **Missing SNP IDs error when using `ld_reference_file`** — LD computation joins on the `rs` column; ensure your GWAS DataFrame has an `rs` (or configured `rs_col`) column populated with SNP IDs matching the PLINK `.bim`.
- **`npx: command not found` on `pre-commit install`** — Contributors only. Install Node.js 20+ (`brew install node` on macOS, `sudo apt install nodejs npm` on Debian/Ubuntu).

## Next Steps

- [ARCHITECTURE.md](ARCHITECTURE.md) — How the plotter, backends, and reference data fit together.
- [CONFIGURATION.md](CONFIGURATION.md) — Plotter options, species config, recombination/LD cache locations.
- [USER_GUIDE.md](USER_GUIDE.md) — Task-oriented walkthroughs (Manhattan, PheWAS, fine-mapping, forest plots).
- [CODEMAP.md](CODEMAP.md) — File-by-file tour of `src/pylocuszoom/`.
- `examples/getting_started.ipynb` — Interactive notebook covering the regional plot workflow.

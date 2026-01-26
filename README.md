# pyLocusZoom

[![CI](https://github.com/michael-denyer/pyLocusZoom/actions/workflows/ci.yml/badge.svg)](https://github.com/michael-denyer/pyLocusZoom/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.5+-11557c.svg)](https://matplotlib.org/)
[![Plotly](https://img.shields.io/badge/Plotly-5.0+-3F4F75.svg)](https://plotly.com/python/)
[![Bokeh](https://img.shields.io/badge/Bokeh-3.8+-E6526F.svg)](https://bokeh.org/)
[![Pandas](https://img.shields.io/badge/Pandas-1.4+-150458.svg)](https://pandas.pydata.org/)

<img src="logo.svg" alt="pyLocusZoom logo" width="120" align="right">

Regional association plots for GWAS results with LD coloring, gene tracks, and recombination rate overlays.

Inspired by [LocusZoom](http://locuszoom.org/) and [locuszoomr](https://github.com/myles-lewis/locuszoomr).

## Features

- **LD coloring**: SNPs colored by linkage disequilibrium (R²) with lead variant
- **Gene track**: Annotated gene/exon positions below the association plot
- **Recombination rate**: Overlay showing recombination rate across region (*Canis lupus familiaris* only)
- **SNP labels**: Automatic labeling of top SNPs with RS ID or nearest gene
- **Species support**: Built-in *Canis lupus familiaris* (CanFam3.1/CanFam4), *Felis catus* (FelCat9), or custom species
- **CanFam4 support**: Automatic coordinate liftover for recombination maps
- **Multiple backends**: matplotlib (static), plotly (interactive), bokeh (dashboards)
- **Stacked plots**: Compare multiple GWAS/phenotypes vertically
- **eQTL overlay**: Expression QTL data as separate panel
- **PySpark support**: Handles large-scale genomics DataFrames

![Example regional association plot](examples/regional_plot.png)

## Installation

```bash
uv add pylocuszoom
```

Or with pip:

```bash
pip install pylocuszoom
```

## Quick Start

```python
from pylocuszoom import LocusZoomPlotter

# Initialize plotter (loads reference data for canine)
plotter = LocusZoomPlotter(species="canine")

# Create regional plot
fig = plotter.plot(
    gwas_df,                    # DataFrame with ps, p_wald, rs columns
    chrom=1,
    start=1000000,
    end=2000000,
    lead_pos=1500000,           # Highlight lead SNP
)

fig.savefig("regional_plot.png", dpi=150)
```

## Full Example

```python
from pylocuszoom import LocusZoomPlotter

plotter = LocusZoomPlotter(
    species="canine",                   # or "feline", or None for custom
    plink_path="/path/to/plink",        # Optional, auto-detects if on PATH
)

fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=1000000,
    end=2000000,
    lead_pos=1500000,
    ld_reference_file="genotypes.bed",  # For LD calculation
    genes_df=genes_df,                  # Gene annotations
    exons_df=exons_df,                  # Exon annotations
    show_recombination=True,            # Overlay recombination rate
    snp_labels=True,                    # Label top SNPs
    label_top_n=5,                      # How many to label
    pos_col="ps",                       # Column name for position
    p_col="p_wald",                     # Column name for p-value
    rs_col="rs",                        # Column name for SNP ID
    figsize=(12, 8),
)
```

## Genome Builds

The default genome build for canine is CanFam3.1. For CanFam4 data:

```python
plotter = LocusZoomPlotter(species="canine", genome_build="canfam4")
```

Recombination maps are automatically lifted over from CanFam3.1 to CanFam4 coordinates using the UCSC liftOver chain file.

## Using with Other Species

```python
# Feline (LD and gene tracks, user provides recombination data)
plotter = LocusZoomPlotter(species="feline")

# Custom species (provide all reference data)
plotter = LocusZoomPlotter(
    species=None,
    recomb_data_dir="/path/to/recomb_maps/",
)

# Or provide data per-plot
fig = plotter.plot(
    gwas_df,
    chrom=1, start=1000000, end=2000000,
    recomb_df=my_recomb_dataframe,
    genes_df=my_genes_df,
)
```

## Interactive Backends (Coming Soon)

> **Note:** Interactive backends (plotly, bokeh) are planned but not yet fully integrated. Currently all plots use matplotlib.

```python
# Static publication-quality plot (default, currently only supported backend)
plotter = LocusZoomPlotter(species="canine", backend="matplotlib")
fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
fig.savefig("plot.png", dpi=150)
```

Future releases will support:
- **Plotly**: Interactive plots with hover tooltips, zoom/pan
- **Bokeh**: Dashboard-friendly interactive plots

## Stacked Plots

Compare multiple GWAS results vertically with shared x-axis:

```python
fig = plotter.plot_stacked(
    [gwas_height, gwas_bmi, gwas_whr],
    chrom=1,
    start=1000000,
    end=2000000,
    panel_labels=["Height", "BMI", "WHR"],
    genes_df=genes_df,
)
```

## eQTL Overlay

Add expression QTL data as a separate panel:

```python
eqtl_df = pd.DataFrame({
    "pos": [1000500, 1001200, 1002000],
    "p_value": [1e-6, 1e-4, 0.01],
    "gene": ["BRCA1", "BRCA1", "BRCA1"],
})

fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1, start=1000000, end=2000000,
    eqtl_df=eqtl_df,
    eqtl_gene="BRCA1",
    genes_df=genes_df,
)
```

## PySpark Support

For large-scale genomics data, pass PySpark DataFrames directly:

```python
from pylocuszoom import LocusZoomPlotter, to_pandas

# PySpark DataFrame (automatically converted)
fig = plotter.plot(spark_gwas_df, chrom=1, start=1000000, end=2000000)

# Or convert manually with sampling for very large data
pandas_df = to_pandas(spark_gwas_df, sample_size=100000)
```

Install PySpark support: `uv add pylocuszoom[spark]`

## Data Formats

### GWAS Results DataFrame

Required columns (names configurable via `pos_col`, `p_col`, `rs_col`):

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `ps` | int | Yes | Genomic position in base pairs (1-based). Must match coordinate system of genes/recombination data. |
| `p_wald` | float | Yes | Association p-value (0 < p ≤ 1). Values are -log10 transformed for plotting. |
| `rs` | str | No | SNP identifier (e.g., "rs12345" or "chr1:12345"). Used for labeling top SNPs if `snp_labels=True`. |

Example:
```python
gwas_df = pd.DataFrame({
    "ps": [1000000, 1000500, 1001000],
    "p_wald": [1e-8, 1e-6, 0.05],
    "rs": ["rs123", "rs456", "rs789"],
})
```

### Genes DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str or int | Yes | Chromosome identifier. Accepts "1", "chr1", or 1. The "chr" prefix is stripped for matching. |
| `start` | int | Yes | Gene start position (bp, 1-based). Transcript start for strand-aware genes. |
| `end` | int | Yes | Gene end position (bp, 1-based). Must be ≥ start. |
| `gene_name` | str | Yes | Gene symbol displayed in track (e.g., "BRCA1", "TP53"). Keep short for readability. |

Example:
```python
genes_df = pd.DataFrame({
    "chr": ["1", "1", "1"],
    "start": [1000000, 1050000, 1100000],
    "end": [1020000, 1080000, 1150000],
    "gene_name": ["GENE1", "GENE2", "GENE3"],
})
```

### Exons DataFrame (optional)

Provides exon/intron structure. If omitted, genes are drawn as simple rectangles.

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str or int | Yes | Chromosome identifier. |
| `start` | int | Yes | Exon start position (bp). |
| `end` | int | Yes | Exon end position (bp). |
| `gene_name` | str | Yes | Parent gene symbol. Must match `gene_name` in genes DataFrame. |

### Recombination DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `pos` | int | Yes | Genomic position (bp). Should span the plotted region with reasonable density (every ~10kb). |
| `rate` | float | Yes | Recombination rate in centiMorgans per megabase (cM/Mb). Typical range: 0-50 cM/Mb. |

Example:
```python
recomb_df = pd.DataFrame({
    "pos": [1000000, 1010000, 1020000],
    "rate": [0.5, 2.3, 1.1],
})
```

### Recombination Map Files

When using `recomb_data_dir`, files must be named `chr{N}_recomb.tsv` (e.g., `chr1_recomb.tsv`, `chrX_recomb.tsv`).

Format: Tab-separated with header row:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome number (without "chr" prefix) |
| `pos` | Position in base pairs |
| `rate` | Recombination rate (cM/Mb) |
| `cM` | Cumulative genetic distance (optional, not used for plotting) |

```
chr	pos	rate	cM
1	10000	0.5	0.005
1	20000	1.2	0.017
1	30000	0.8	0.025
```

## Reference Data

Canine recombination maps are downloaded from [Campbell et al. 2016](https://github.com/cflerin/dog_recombination) on first use.

To manually download:

```python
from pylocuszoom import download_canine_recombination_maps

download_canine_recombination_maps()
```

## Logging

Logging uses [loguru](https://github.com/Delgan/loguru) and is configured via the `log_level` parameter (default: `"INFO"`):

```python
# Suppress logging
plotter = LocusZoomPlotter(log_level=None)

# Enable DEBUG level for troubleshooting
plotter = LocusZoomPlotter(log_level="DEBUG")
```

## Requirements

- Python >= 3.10
- matplotlib >= 3.5.0
- pandas >= 1.4.0
- numpy >= 1.21.0
- loguru >= 0.7.0
- plotly >= 5.0.0
- bokeh >= 3.8.2
- kaleido >= 0.2.0 (for plotly static export)
- pyliftover >= 0.4 (for CanFam4 coordinate liftover)
- [PLINK 1.9](https://www.cog-genomics.org/plink/) (for LD calculations) - must be on PATH or specify `plink_path`

Optional:
- pyspark >= 3.0.0 (for PySpark DataFrame support) - `uv add pylocuszoom[spark]`

## License

GPL-3.0-or-later

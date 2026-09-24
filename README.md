[![CI](https://github.com/michael-denyer/pyLocusZoom/actions/workflows/ci.yml/badge.svg)](https://github.com/michael-denyer/pyLocusZoom/actions/workflows/ci.yml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22665975-orange)](https://doi.org/10.5281/zenodo.22665975)
[![PyPI](https://img.shields.io/pypi/v/pylocuszoom)](https://pypi.org/project/pylocuszoom/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-red.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.5+-11557c.svg)](https://matplotlib.org/)
[![Plotly](https://img.shields.io/badge/Plotly-5.15+-3F4F75.svg)](https://plotly.com/python/)
[![Bokeh](https://img.shields.io/badge/Bokeh-3.8+-E6526F.svg)](https://bokeh.org/)
[![Buy Me A Coffee](https://img.shields.io/badge/Buy%20Me%20A%20Coffee-support-yellow?logo=buy-me-a-coffee&logoColor=white)](https://buymeacoffee.com/codenyer)
<img src="https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/logo.svg" alt="pyLocusZoom logo" width="120" align="right">
# pyLocusZoom

Designed for publication-ready GWAS visualization with regional association plots, gene tracks, eQTL, PheWAS, fine-mapping, and forest plots.

Inspired by [LocusZoom](http://locuszoom.org/) and [locuszoomr](https://github.com/myles-lewis/locuszoomr).

## Features

1. **Regional association plot**:

    - **Multi-species support**: Built-in reference data for *Canis lupus familiaris* (CanFam3.1/CanFam4) and *Felis catus* (FelCat9), or optionally provide your own for any species
    - **LD coloring**: SNPs colored by linkage disequilibrium (R²) with lead variant
    - **Gene tracks**: Annotated gene/exon positions below the association plot
    - **Recombination rate**: Overlay across region (*Canis lupus familiaris* built-in, or user-provided)
    - **SNP labels (matplotlib)**: Automatic labeling of top SNPs by p-value (RS IDs)
    - **Hover tooltips (Plotly and Bokeh)**: Detailed SNP data on hover

    ![Example regional association plot with LD coloring, gene track, and recombination overlay](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/regional_plot_with_recomb.png)
    *Regional association plot with LD coloring, gene/exon track, recombination rate overlay (blue line), and top SNP labels.*

2. **Stacked plots**: Compare multiple GWAS/phenotypes vertically
3. **Miami plots**: Mirrored Manhattan plots for comparing two GWAS datasets (discovery vs replication)
4. **Manhattan plots**: Genome-wide association visualization with chromosome coloring
5. **QQ plots**: Quantile-quantile plots with confidence bands and genomic inflation factor
6. **eQTL plot**: Expression QTL data aligned with association plots and gene tracks
7. **Fine-mapping plots**: Visualize SuSiE credible sets with posterior inclusion probabilities
8. **PheWAS plots**: Phenome-wide association study visualization across multiple phenotypes
9. **Forest plots**: Meta-analysis effect size visualization with confidence intervals
10. **LD heatmaps**: Triangular heatmaps showing pairwise LD patterns, standalone or integrated below regional plots
11. **Colocalization plots**: GWAS-eQTL scatter plots with LD coloring, correlation statistics, and effect direction visualization
12. **Multiple backends**: matplotlib (publication-ready), plotly (interactive), bokeh (dashboard integration)
13. **Pandas and PySpark support**: Works with both Pandas and PySpark DataFrames for large-scale genomics data
14. **Convenience data file loaders**: Load and validate common GWAS, eQTL and fine-mapping file formats
15. **Automatic gene annotations**: Fetch gene/exon data with caching, from UCSC for CanFam3.1, CanFam4 and FelCat9 and from the Ensembl REST API for human, mouse, rat and any other Ensembl species

## Installation

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

PySpark DataFrame support is an extra: `pip install "pylocuszoom[spark]"`.
LD colouring from a genotype fileset needs [PLINK 1.9](https://www.cog-genomics.org/plink/)
on your PATH, or its location passed as `plink_path`.

## Quick Start

```python
from pylocuszoom import LDConfig, LocusZoomPlotter

# Initialize plotter (loads reference data for canine)
plotter = LocusZoomPlotter(species="canine", auto_genes=True)

# The region is passed directly; every other option lives on a config model
fig = plotter.plot(
    gwas_df,                        # DataFrame with chr, pos, p_value, rs columns
    chrom=1,
    start=1000000,
    end=2000000,
    ld=LDConfig(lead_pos=1500000),  # Highlight lead SNP
)
fig.savefig("regional_plot.png", dpi=150)
```

Pass `backend="plotly"` or `backend="bokeh"` to the plotter for an interactive
figure with hover tooltips. The [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md) covers every plot type, the
config models, the file loaders, the input column formats and species support.

## Gallery

Each figure below links to the User Guide section with its code.

### Stacked regional plots

Compare several GWAS over one region, with a shared gene track. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#stacked-plots)

![Example stacked plot comparing two phenotypes](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/stacked_plot.png)
*Stacked plot comparing two phenotypes with LD coloring and shared gene track.*

### eQTL overlay

Expression QTL results in their own panel below the association plot. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#eqtl-overlay)

![Example eQTL overlay plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/eqtl_overlay.png)
*eQTL overlay with effect direction (up/down triangles) and magnitude binning.*

### Fine-mapping

SuSiE, FINEMAP, CAVIAR or PolyFun results, with credible sets coloured. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#fine-mapping-visualization)

![Example fine-mapping plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/finemapping_plot.png)
*Fine-mapping visualization with PIP line and credible set coloring (CS1/CS2).*

### LD heatmaps

Triangular pairwise LD heatmaps, standalone or as a panel below a regional plot. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#ld-heatmaps)

![Example LD heatmap](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/ld_heatmap.png)
*Triangular LD heatmap with R² values and lead SNP highlighted.*

![Example regional plot with LD heatmap](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/regional_with_ld_heatmap.png)
*Regional association plot with integrated LD heatmap panel below.*

### Colocalization

GWAS against eQTL significance, coloured by LD or by effect-direction agreement. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#colocalization-plots)

![Example colocalization plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/colocalization_plot.png)
*GWAS-eQTL colocalization scatter plot with LD coloring and correlation statistics.*

### PheWAS

One variant's associations across phenotypes, grouped by category. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#phewas-plots)

![Example PheWAS plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/phewas_plot.png)
*PheWAS plot showing associations across phenotype categories with significance threshold.*

### Forest plots

Effect sizes with confidence intervals across studies. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#forest-plots)

![Example forest plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/forest_plot.png)
*Forest plot with effect sizes, confidence intervals, and weight-proportional markers.*

### Miami plots

Two GWAS mirrored about a shared chromosome axis. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#miami-plots)

![Example Miami plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/miami_plot.png)
*Miami plot comparing discovery and replication GWAS with mirrored y-axes and region highlighting.*

### Manhattan plots

Genome-wide associations by chromosome, or by category. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#manhattan-plots)

![Example Manhattan plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/manhattan_plot.png)
*Manhattan plot showing genome-wide associations with chromosome coloring and significance threshold.*

### QQ plots

Observed against expected p-values, with a confidence band and λ. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#qq-plots)

![Example QQ plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/qq_plot.png)
*QQ plot with 95% confidence band and genomic inflation factor (λ).*

### Stacked Manhattan plots

Several GWAS on one chromosome axis. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#stacked-manhattan-plots)

![Example stacked Manhattan plot](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/manhattan_stacked.png)
*Stacked Manhattan plots comparing three GWAS studies with shared chromosome axis.*

### Manhattan and QQ side by side

A one-figure GWAS summary. Every genome-wide plot takes a `GenomeWideStyle` for palette, point and font styling. [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md#manhattan-and-qq-side-by-side)

![Example Manhattan and QQ side-by-side](https://raw.githubusercontent.com/michael-denyer/pyLocusZoom/main/examples/matplotlib/manhattan_qq_sidebyside.png)
*Combined Manhattan and QQ plot showing genome-wide associations and p-value distribution.*

## Documentation

- [Getting Started](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/GETTING-STARTED.md) - Installation and first plot
- [User Guide](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/USER_GUIDE.md) - Every plot type, the config models, loaders, data formats and species support
- [Migrating to 5.0](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/MIGRATING-5.0.md) - What to change when upgrading from 4.x
- [Configuration](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/CONFIGURATION.md) - Cache locations and the environment variables that move them
- [Architecture](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/ARCHITECTURE.md) - Design decisions and component overview
- [Code Map](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/CODEMAP.md) - Architecture diagram with source code links
- [Development](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/DEVELOPMENT.md) - Dev setup, pre-commit hooks, contributing workflow
- [Testing](https://github.com/michael-denyer/pyLocusZoom/blob/main/docs/TESTING.md) - Running and writing tests
- [Example Notebook](https://github.com/michael-denyer/pyLocusZoom/blob/main/examples/getting_started.ipynb) - Interactive tutorial
- [CHANGELOG](https://github.com/michael-denyer/pyLocusZoom/blob/main/CHANGELOG.md) - Version history

## Citation

If you use pyLocusZoom in your research, please cite it. GitHub's "Cite this repository" button reads [CITATION.cff](https://github.com/michael-denyer/pyLocusZoom/blob/main/CITATION.cff), and each GitHub release is archived on Zenodo with its own DOI. The concept DOI [10.5281/zenodo.22665975](https://doi.org/10.5281/zenodo.22665975) always resolves to the latest version.

```bibtex
@software{denyer_pylocuszoom,
  author  = {Denyer, Michael},
  title   = {pyLocusZoom: Python library for multi-species GWAS visualization},
  url     = {https://github.com/michael-denyer/pyLocusZoom},
  doi     = {10.5281/zenodo.22665975},
  license = {GPL-3.0-or-later}
}
```

## License

GPL-3.0-or-later

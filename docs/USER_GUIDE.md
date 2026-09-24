# pyLocusZoom User Guide

[![PyPI version](https://img.shields.io/pypi/v/pylocuszoom.svg)](https://pypi.org/project/pylocuszoom/)

Comprehensive documentation for pyLocusZoom - regional association plots for GWAS results with LD coloring, gene tracks, and recombination rate overlays.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Plot Types](#plot-types)
  - [Regional Association Plot](#regional-association-plot)
  - [Stacked Plots](#stacked-plots)
  - [eQTL Overlay](#eqtl-overlay)
  - [Fine-mapping Visualization](#fine-mapping-visualization)
  - [LD Heatmaps](#ld-heatmaps)
  - [Colocalization Plots](#colocalization-plots)
  - [PheWAS Plots](#phewas-plots)
  - [Forest Plots](#forest-plots)
  - [Miami Plots](#miami-plots)
  - [Manhattan Plots](#manhattan-plots)
  - [QQ Plots](#qq-plots)
  - [Stacked Manhattan Plots](#stacked-manhattan-plots)
  - [Manhattan and QQ Side-by-Side](#manhattan-and-qq-side-by-side)
- [Backends](#backends)
  - [Matplotlib (Static)](#matplotlib-static)
  - [Plotly (Interactive)](#plotly-interactive)
  - [Bokeh (Dashboard)](#bokeh-dashboard)
  - [Custom Backends](#custom-backends)
- [Plotter Reference](#plotter-reference)
  - [LocusZoomPlotter](#locuszoomplotter)
  - [plot() Method](#plot-method)
  - [plot_stacked() Method](#plot_stacked-method)
- [File Loaders](#file-loaders)
  - [GWAS Loaders](#gwas-loaders)
  - [eQTL Loaders](#eqtl-loaders)
  - [Fine-mapping Loaders](#fine-mapping-loaders)
  - [Gene Annotation Loaders](#gene-annotation-loaders)
- [Data Formats](#data-formats)
- [Species Support](#species-support)
- [Recipes & Examples](#recipes--examples)
- [API Stability](#api-stability)

---

## Installation

### pip (PyPI)

```bash
pip install pylocuszoom
```

### uv

```bash
uv add pylocuszoom
```

### conda (Bioconda)

```bash
conda install -c bioconda pylocuszoom
```

### Optional Dependencies

```bash
# For PySpark DataFrame support
pip install pylocuszoom[spark]
```

### External Requirements

**PLINK 1.9** is required for LD calculations. Install from [cog-genomics.org/plink](https://www.cog-genomics.org/plink/) and ensure it's on your PATH, or specify the path via `plink_path` parameter.

---

## Quick Start

```python
import pandas as pd
from pylocuszoom import LDConfig, LocusZoomPlotter

# Sample GWAS data
gwas_df = pd.DataFrame({
    "chr": [1] * 5,
    "pos": [1000000, 1000500, 1001000, 1001500, 1002000],
    "p_value": [0.05, 1e-4, 1e-8, 1e-6, 0.01],
    "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
})

# Create plotter
plotter = LocusZoomPlotter(species="canine")

# The region is passed directly; every other option lives on a config model
fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=999000,
    end=1003000,
    ld=LDConfig(lead_pos=1001000),  # Highlight the most significant SNP
)

# Save
fig.savefig("my_plot.png", dpi=150, bbox_inches="tight")
```

---

## Plot Types

### Regional Association Plot

The basic single-panel plot showing association signals with LD coloring.

![Regional association plot](../examples/matplotlib/regional_plot.png)

```python
from pylocuszoom import DisplayConfig, LDConfig, PanelInputs

fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=1000000,
    end=2000000,
    ld=LDConfig(
        lead_pos=1500000,           # Lead SNP position
        ld_reference_file="geno",   # PLINK fileset for LD calculation
    ),
    display=DisplayConfig(
        show_recombination=True,    # Recombination rate overlay
        snp_labels=True,            # Label top SNPs (matplotlib only)
        label_top_n=5,              # How many to label
    ),
    panels=PanelInputs(
        genes_df=genes_df,          # Gene annotations
        exons_df=exons_df,          # Exon structure (optional)
    ),
)
```

**Features:**

- SNPs coloured by R² with the lead variant in five bins: blue below 0.2, then
  cyan, green, orange, and red from 0.8; grey where R² is unknown
- Lead SNP shown as purple diamond
- Gene track with intron/exon structure
- Recombination rate overlay (blue line, right y-axis)
- Genome-wide significance line (red dashed, default 5e-8)
- Top SNPs labelled with their ids from the `rs_col` column (matplotlib only)

### Stacked Plots

Compare multiple GWAS results vertically with a shared x-axis.

![Stacked plot](../examples/matplotlib/stacked_plot.png)

```python
fig = plotter.plot_stacked(
    [gwas_height, gwas_bmi, gwas_whr],
    chrom=1,
    start=1000000,
    end=2000000,
    lead_positions=[1500000, 1480000, 1520000],  # Per-panel leads
    panel_labels=["Height", "BMI", "WHR"],
    panels=PanelInputs(genes_df=genes_df),
)
```

**Features:**

- Vertical stacking with aligned x-axes
- Independent LD coloring per panel
- Shared gene track at bottom
- Optional recombination overlay (top panel only)

### eQTL Overlay

Add expression QTL data as a separate panel below GWAS results.

![eQTL overlay](../examples/matplotlib/eqtl_overlay.png)

```python
from pylocuszoom import EqtlInput, PanelInputs

eqtl_df = pd.DataFrame(
    {
        "chr": [1] * 3,
        "pos": [1000500, 1001200, 1002000],
        "p_value": [1e-6, 1e-4, 0.01],
        "gene": ["BRCA1", "BRCA1", "BRCA1"],
        "effect_size": [0.5, -0.3, 0.1],  # Optional: colors by effect direction
    }
)

fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1000000,
    end=2000000,
    panels=PanelInputs(genes_df=genes_df, eqtl=EqtlInput(data=eqtl_df, gene="BRCA1")),
)
```

**Features:**

- Separate panel for eQTL associations
- Color by effect direction and size (warm up-triangles = positive, cool
  down-triangles = negative)
- Filter to specific target gene

### Fine-mapping Visualization

Visualize SuSiE or other fine-mapping results with credible set coloring.

![Fine-mapping plot](../examples/matplotlib/finemapping_plot.png)

```python
from pylocuszoom import FinemappingInput, PanelInputs

finemapping_df = pd.DataFrame(
    {
        "chr": [1] * 4,
        "pos": [1000500, 1001200, 1002000, 1003500],
        "pip": [0.85, 0.12, 0.02, 0.45],  # Posterior inclusion probability
        "cs": [1, 1, 0, 2],  # Credible set (0 = not in CS)
    }
)

fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1000000,
    end=2000000,
    panels=PanelInputs(
        genes_df=genes_df,
        finemapping=FinemappingInput(data=finemapping_df, cs_col="cs"),
    ),
)
```

**Features:**

- PIP values shown as line plot
- Credible sets colored distinctly (CS1 = orange, CS2 = blue, CS3 = green, etc.)
- Variants not in credible sets shown in gray

### LD Heatmaps

Create triangular LD heatmaps showing pairwise linkage disequilibrium patterns.

![LD heatmap](../examples/matplotlib/ld_heatmap.png)

```python
from pylocuszoom import LDHeatmapPlotter

# ld_matrix is a square DataFrame with SNP IDs as index/columns
# snp_ids is a list of SNP IDs in matrix order

ld_plotter = LDHeatmapPlotter()
fig = ld_plotter.plot_ld_heatmap(
    ld_matrix,
    snp_ids=snp_ids,
    lead_snp="rs12345",    # Highlight lead SNP (red)
    highlight_snps=None,   # Optional extra SNPs to mark (blue)
    metric="r2",           # or "dprime"
)
fig.savefig("ld_heatmap.png", dpi=150)
```

**Features:**

- Triangular heatmap showing pairwise LD (R² or D')
- White-to-red color gradient
- Optional lead SNP highlighting
- Colorbar legend with metric label

#### Integrated LD Heatmap with Regional Plot

Add an LD heatmap panel below a regional association plot:

![Regional plot with LD heatmap](../examples/matplotlib/regional_with_ld_heatmap.png)

```python
from pylocuszoom import LDConfig, LDHeatmapInput, LocusZoomPlotter, PanelInputs

plotter = LocusZoomPlotter(species="canine")

fig = plotter.plot(
    gwas_df,
    chrom=1,
    start=1000000,
    end=2000000,
    ld=LDConfig(lead_pos=1500000),
    panels=PanelInputs(
        ld_heatmap=LDHeatmapInput(matrix=ld_matrix, snp_ids=snp_ids, height=0.25)
    ),
)
```

**Features:**

- Heatmap panel automatically added below association plot
- SNPs aligned with x-axis coordinates from GWAS data
- Works with both `plot()` and `plot_stacked()`
- Heatmap at very bottom in stacked plots

### Colocalization Plots

Visualize GWAS-eQTL colocalization by comparing association signals in a scatter plot with LD coloring.

![Colocalization plot](../examples/matplotlib/colocalization_plot.png)

```python
from pylocuszoom import ColocConfig, ColocPlotter
import pandas as pd

# GWAS data with position and p-value
gwas_df = pd.DataFrame({
    "pos": [1000000, 1000500, 1001000, 1001500, 1002000],
    "p": [1e-8, 1e-6, 1e-4, 0.01, 0.05],
    "ld_r2": [1.0, 0.8, 0.5, 0.2, 0.1],  # Optional: LD with lead SNP
})

# eQTL data with position and p-value
eqtl_df = pd.DataFrame({
    "pos": [1000000, 1000500, 1001000, 1001500, 1002000],
    "p": [1e-6, 1e-8, 1e-5, 0.02, 0.1],
})

plotter = ColocPlotter()
fig = plotter.plot_coloc(
    gwas_df,
    eqtl_df,
    # ld_col is optional: it colours the points by LD
    config=ColocConfig(pos_col="pos", gwas_p_col="p", eqtl_p_col="p", ld_col="ld_r2"),
)
fig.savefig("colocalization.png", dpi=150)
```

**Features:**

- Scatter plot comparing GWAS -log10(p) vs eQTL -log10(p)
- Points colored by LD (R²) with the lead SNP, in the regional plot's five bins
- Lead SNP labeled on plot
- Pearson correlation coefficient and p-value displayed
- Significance threshold reference lines

**Effect Direction Coloring:**

Color points by whether GWAS and eQTL effects are in the same direction (congruent) or opposite (incongruent). This helps identify whether increased gene expression is associated with increased or decreased disease risk.

![Colocalization effect plot](../examples/matplotlib/colocalization_effect_plot.png)

```python
# Effect direction coloring requires effect size columns in both datasets
gwas_df = pd.DataFrame({
    "pos": [1000000, 1000500, 1001000],
    "p": [1e-8, 1e-6, 1e-4],
    "beta": [0.5, 0.3, -0.2],  # GWAS effect sizes
})

eqtl_df = pd.DataFrame({
    "pos": [1000000, 1000500, 1001000],
    "p": [1e-6, 1e-8, 1e-5],
    "slope": [0.8, 0.5, -0.4],  # eQTL effect sizes
})

fig = plotter.plot_coloc(
    gwas_df,
    eqtl_df,
    config=ColocConfig(
        pos_col="pos",
        gwas_p_col="p",
        eqtl_p_col="p",
        gwas_effect_col="beta",
        eqtl_effect_col="slope",
        color_by_effect=True,  # Green=same direction, Red=opposite
        h4_posterior=0.85,  # Display coloc H4 posterior probability
    ),
)
```

**Additional options:**

```python
fig = plotter.plot_coloc(
    gwas_df,
    eqtl_df,
    config=ColocConfig(pos_col="pos", gwas_p_col="p", eqtl_p_col="p", ld_col="ld_r2"),
    # Significance thresholds, per call; None draws no line
    gwas_threshold=5e-8,
    eqtl_threshold=1e-5,
)
```

**Interactive Plotly backend:**

```python
from pylocuszoom import ColocConfig, ColocPlotter

plotter = ColocPlotter(backend="plotly")
fig = plotter.plot_coloc(
    gwas_df,
    eqtl_df,
    config=ColocConfig(pos_col="pos", gwas_p_col="p", eqtl_p_col="p", ld_col="ld_r2"),
)

# Save as interactive HTML
fig.write_html("colocalization_interactive.html")

# Display in Jupyter notebook
from IPython.display import display, HTML
display(HTML(fig.to_html(include_plotlyjs='cdn')))
```

### PheWAS Plots

Visualize associations of a single variant across multiple phenotypes in a phenome-wide association study.

![PheWAS plot](../examples/matplotlib/phewas_plot.png)

```python
from pylocuszoom import StatsPlotter

phewas_df = pd.DataFrame({
    "phenotype": ["Height", "BMI", "T2D", "CAD", "HDL"],
    "p_value": [1e-15, 0.05, 1e-8, 1e-3, 1e-10],
    "category": ["Anthropometric", "Anthropometric", "Metabolic", "Cardiovascular", "Lipids"],
})

plotter = StatsPlotter()
fig = plotter.plot_phewas(
    phewas_df,
    variant_id="rs12345",
    category_col="category",
    significance_threshold=5e-8,
)
```

**Features:**

- Phenotypes grouped and colored by category
- Genome-wide significance line (red dashed)
- Optional effect direction markers (triangles for +/-)
- 12-color palette for distinct categories

### Forest Plots

Create forest plots for meta-analysis visualization showing effect sizes with confidence intervals.

![Forest plot](../examples/matplotlib/forest_plot.png)

```python
from pylocuszoom import StatsPlotter

forest_df = pd.DataFrame({
    "study": ["Study A", "Study B", "Study C", "Meta-analysis"],
    "effect": [0.45, 0.52, 0.38, 0.46],
    "ci_lower": [0.30, 0.35, 0.20, 0.40],
    "ci_upper": [0.60, 0.69, 0.56, 0.52],
    "weight": [25, 35, 20, 100],  # Optional: affects marker size
})

plotter = StatsPlotter()
fig = plotter.plot_forest(
    forest_df,
    variant_id="rs12345",
    weight_col="weight",
    null_value=0.0,  # Reference line (0 for beta, 1 for OR)
    effect_label="Effect Size",
)
```

**Features:**

- Effect sizes as squares with confidence interval whiskers
- Marker size scaled by study weight (optional)
- Null effect reference line
- Study names as y-axis labels

### Miami Plots

Miami plots (mirrored Manhattan plots) compare two GWAS datasets side-by-side with a shared x-axis. The top panel shows -log10(p) ascending, while the bottom panel is inverted.

![Miami plot](../examples/matplotlib/miami_plot.png)

```python
import pandas as pd

from pylocuszoom import MiamiPlotter

# Two GWAS datasets to compare
gwas1 = pd.read_csv("gwas_study1.csv")
gwas2 = pd.read_csv("gwas_study2.csv")

plotter = MiamiPlotter()
fig = plotter.plot_miami(
    top_df=gwas1,
    bottom_df=gwas2,
    top_label="Study 1",
    bottom_label="Study 2",
    figsize=(14, 8),
)
fig.savefig("miami.png", dpi=150)
```

**Features:**

- Mirrored panels with shared x-axis and consistent chromosome colors
- Per-panel significance thresholds (`top_threshold`, `bottom_threshold`)
- Panel labels to identify datasets
- SNP annotations independent per panel (`top_snp_annotations`, `bottom_snp_annotations`)
- Region highlighting across both panels (`highlight_regions`)
- Interactive hover tooltips in plotly/bokeh backends
- Full support for all three backends (matplotlib, plotly, bokeh)

**Customization options:**

```python
fig = plotter.plot_miami(
    top_df=gwas1,
    bottom_df=gwas2,
    # Per-panel thresholds
    top_threshold=5e-8,
    bottom_threshold=1e-5,
    # Panel labels
    top_label="Discovery Cohort",
    bottom_label="Replication Cohort",
    # SNP annotations (list of SNP IDs; requires rs_col to be set)
    rs_col="rs",
    top_snp_annotations=["rs123"],
    bottom_snp_annotations=["rs456"],
    # Highlight regions (list of (chrom, start, end) tuples)
    highlight_regions=[(1, 1000000, 2000000), (5, 50000000, 51000000)],
)
```

**Interactive Plotly backend:**

```python
from pylocuszoom import MiamiPlotter

# Create plotter with Plotly backend
plotter = MiamiPlotter(species="human", backend="plotly")

fig = plotter.plot_miami(
    top_df=discovery_df,
    bottom_df=replication_df,
    top_label="Discovery",
    bottom_label="Replication",
    top_threshold=5e-8,
    bottom_threshold=5e-8,
    highlight_regions=[("6", 25_000_000, 35_000_000)],  # MHC region
    title="Discovery vs Replication GWAS",
)

# Save as interactive HTML
fig.write_html("miami_interactive.html")

# Display in Jupyter notebook
from IPython.display import display, HTML
display(HTML(fig.to_html(include_plotlyjs='cdn')))
```

**Interactive Bokeh backend:**

```python
from pylocuszoom import MiamiPlotter
from bokeh.io import output_file, save
from bokeh.resources import CDN
from bokeh.embed import file_html

# Create plotter with Bokeh backend
plotter = MiamiPlotter(species="human", backend="bokeh")

fig = plotter.plot_miami(
    top_df=discovery_df,
    bottom_df=replication_df,
    top_label="Discovery",
    bottom_label="Replication",
    top_threshold=5e-8,
    bottom_threshold=5e-8,
    highlight_regions=[("6", 25_000_000, 35_000_000)],
    title="Discovery vs Replication GWAS",
)

# Save as interactive HTML
output_file("miami_bokeh.html")
save(fig)

# Display in Jupyter notebook
from IPython.display import display, HTML
display(HTML(file_html(fig, CDN, "Miami Plot")))
```

### Manhattan Plots

Genome-wide Manhattan plots showing associations across all chromosomes.

![Manhattan plot](../examples/matplotlib/manhattan_plot.png)

```python
from pylocuszoom import ManhattanPlotter

plotter = ManhattanPlotter()
fig = plotter.plot_manhattan(
    gwas_df,
    significance_threshold=5e-8,
    figsize=(12, 5),
)
fig.savefig("manhattan.png", dpi=150)
```

**Features:**

- Chromosomes colored alternately for distinction
- Genome-wide significance threshold line (red dashed)
- Automatic cumulative position calculation
- Chromosome labels on x-axis

**Categorical Manhattan plots** (PheWAS-style) put a category on the x axis
instead of genomic position:

```python
from pylocuszoom import GenomeWideConfig, ManhattanPlotter

phewas_df = pd.DataFrame({
    "phenotype": ["Height", "BMI", "T2D", "CAD", "HDL"],
    "pvalue": [1e-15, 0.05, 1e-8, 1e-3, 1e-10],
    "phenotype_category": ["Anthropometric", "Anthropometric", "Metabolic", "Cardiovascular", "Lipids"],
})

fig = ManhattanPlotter().plot_manhattan(
    phewas_df,
    category_col="phenotype_category",
    config=GenomeWideConfig(p_col="pvalue"),
)
```

### QQ Plots

Quantile-quantile plots for assessing p-value distribution and detecting systematic bias.

![QQ plot](../examples/matplotlib/qq_plot.png)

```python
from pylocuszoom import ManhattanPlotter

plotter = ManhattanPlotter()
fig = plotter.plot_qq(
    gwas_df,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(6, 6),
)
fig.savefig("qq_plot.png", dpi=150)
```

**Features:**

- Expected vs observed -log10(p) values
- 95% confidence band (beta distribution)
- Genomic inflation factor (λ) in title
- Identity line for reference

### Stacked Manhattan Plots

Compare multiple GWAS studies in vertically stacked Manhattan plots with shared chromosome axis.

![Stacked Manhattan plot](../examples/matplotlib/manhattan_stacked.png)

```python
from pylocuszoom import ManhattanPlotter

plotter = ManhattanPlotter()
fig = plotter.plot_manhattan_stacked(
    [gwas_study1, gwas_study2, gwas_study3],
    panel_labels=["Study 1", "Study 2", "Study 3"],
    significance_threshold=5e-8,
    figsize=(12, 8),
    title="Multi-study GWAS Comparison",
)
fig.savefig("manhattan_stacked.png", dpi=150)
```

**Features:**

- Vertically stacked panels with aligned x-axes
- Shared chromosome coloring across panels
- Independent y-axes per panel
- Panel labels for study identification
- Optional overall figure title

### Manhattan and QQ Side-by-Side

Combined Manhattan and QQ plots in a single figure for comprehensive GWAS summary.

![Manhattan and QQ side-by-side](../examples/matplotlib/manhattan_qq_sidebyside.png)

```python
from pylocuszoom import ManhattanPlotter

plotter = ManhattanPlotter()
fig = plotter.plot_manhattan_qq(
    gwas_df,
    significance_threshold=5e-8,
    show_confidence_band=True,
    show_lambda=True,
    figsize=(14, 5),
    title="GWAS Results",
)
fig.savefig("manhattan_qq.png", dpi=150)
```

**Features:**

- Manhattan plot on left (wider), QQ plot on right
- Shared significance threshold line on Manhattan
- Confidence band and λ on QQ plot
- Optional overall figure title

Three optional arguments, all off by default:

| Argument | Default | Effect |
|----------|---------|--------|
| `suggestive_threshold` | None | P-value for a second, blue dashed line on the Manhattan panel, such as `1e-5`. |
| `lambda_gc` | None | Inflation factor shown on the QQ panel instead of the one computed from the plotted p-values. |
| `footer` | None | One small italic grey line under both panels, such as the sample size or model. |

```python
fig = plotter.plot_manhattan_qq(
    gwas_df,
    significance_threshold=0.05 / len(gwas_df),
    suggestive_threshold=1e-5,
    lambda_gc=lambda_from_pipeline,
    footer="n = 2,345 dogs | GEMMA LMM",
)
```

### Styling genome-wide plots

Every Manhattan, QQ, stacked, side-by-side and Miami method takes
`style=GenomeWideStyle(...)` for the chromosome palette, points, font sizes
and chromosome axis. Fields you leave unset keep the method's current look;
see [GenomeWideStyle](#genomewidestyle) for each field.

```python
import colorcet as cc
from pylocuszoom import GenomeWideStyle, ManhattanPlotter

style = GenomeWideStyle(
    palette=cc.glasbey_dark[::6],  # any list of colours; for a colormap, pass cmap.colors
    point_size=75,
    title_fontsize=30,
    panel_title_fontsize=26,
    axis_label_fontsize=20,
    tick_label_fontsize=18,
    tick_step=2,         # label chromosomes 1, 3, 5, ...
    tick_rotation=90,
)
fig = ManhattanPlotter(species="canine").plot_manhattan_qq(
    gwas_df, title="Coat colour", suggestive_threshold=1e-5, style=style
)
```

---

## Backends

pyLocusZoom supports three rendering backends for different use cases.

| Backend | Output | Best For | SNP Labels |
|---------|--------|----------|------------|
| `matplotlib` | PNG, PDF, SVG | Publications, presentations | Yes (adjustText) |
| `plotly` | Interactive HTML | Web reports, exploration | No (hover instead) |
| `bokeh` | Interactive HTML | Dashboards, web apps | No (hover instead) |

### Matplotlib (Static)

Default backend for publication-quality static plots.

```python
plotter = LocusZoomPlotter(species="canine", backend="matplotlib")
fig = plotter.plot(gwas_df, chrom=1, start=1e6, end=2e6)
fig.savefig("plot.png", dpi=300, bbox_inches="tight")
fig.savefig("plot.pdf")  # Vector format for publications
```

**Unique features:**

- SNP labels with automatic positioning (adjustText library)
- High DPI for print quality
- Vector formats (PDF, SVG) supported

### Plotly (Interactive)

Interactive plots for web reports and data exploration.

```python
plotter = LocusZoomPlotter(species="canine", backend="plotly")
fig = plotter.plot(gwas_df, chrom=1, start=1e6, end=2e6)
fig.write_html("plot.html")
fig.show()  # Opens in browser
```

**Unique features:**

- Hover tooltips showing SNP ID, position, p-value, LD
- Pan and zoom
- Export to PNG/SVG from browser

### Bokeh (Dashboard)

Interactive plots optimized for dashboard integration.

```python
from bokeh.io import output_file, save
from pylocuszoom import LocusZoomPlotter

plotter = LocusZoomPlotter(species="canine", backend="bokeh")
fig = plotter.plot(gwas_df, chrom=1, start=1e6, end=2e6)
output_file("plot.html")
save(fig)
```

**Unique features:**

- Hover tooltips
- Pan and zoom
- Easy integration with Bokeh server applications

### Custom Backends

A custom backend implements the `PlotBackend` protocol in
`pylocuszoom/backends/base.py`, which carries drawing primitives only: legends
and the recombination overlay are composed above it in
`backends/composition.py`. `SupportsSNPLabels` (matplotlib-style repositioned
labels) is the one optional capability, negotiated with a `@runtime_checkable`
protocol: a backend opts in by implementing `add_snp_labels` and out by
omitting it, and still renders every plot family without it. 5.0 changed the
protocol; [MIGRATING-5.0.md](MIGRATING-5.0.md#custom-backends) lists each
signature change, and [ARCHITECTURE.md](ARCHITECTURE.md) describes the seam.

---

## Plotter Reference

This guide is curated. It covers the plotter classes and the parameters most
people reach for, not every public name. For the complete public surface see
the Public API Surface table in [CODEMAP.md](CODEMAP.md).

### Specialized Plotter Classes

pyLocusZoom provides specialized plotter classes for different plot types:

| Class | Purpose |
|-------|---------|
| `LocusZoomPlotter` | Regional association plots with LD coloring |
| `ManhattanPlotter` | Genome-wide Manhattan and QQ plots |
| `MiamiPlotter` | Two-trait mirrored Manhattan plots |
| `StatsPlotter` | PheWAS and forest plots |
| `LDHeatmapPlotter` | Pairwise LD heatmaps |
| `ColocPlotter` | GWAS × eQTL colocalisation scatter |

```python
from pylocuszoom import LocusZoomPlotter, ManhattanPlotter, StatsPlotter

# Regional plots
regional = LocusZoomPlotter(species="canine")
fig = regional.plot(gwas_df, chrom=1, start=1e6, end=2e6)

# Manhattan/QQ plots
manhattan = ManhattanPlotter()
fig = manhattan.plot_manhattan(gwas_df)

# PheWAS/forest plots
stats = StatsPlotter()
fig = stats.plot_phewas(phewas_df, variant_id="rs12345")
```

### LocusZoomPlotter

The main class for creating regional association plots.

```python
plotter = LocusZoomPlotter(
    species="canine",           # "canine", "feline", "human", etc.
    genome_build="canfam3.1",   # Build for coordinate system (auto-selected)
    backend="matplotlib",       # "matplotlib", "plotly", "bokeh"
    plink_path=None,            # Path to PLINK (auto-detects)
    recomb_data_dir=None,       # Custom recombination maps
    genomewide_threshold=5e-8,  # Significance line threshold
    auto_genes=False,           # Auto-fetch gene track from Ensembl
)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `species` | str | `"canine"` | Species name. `"canine"` has built-in recombination maps. |
| `genome_build` | str | Auto | Genome build (`"canfam3.1"`, `"canfam4"`, `"felCat9"`). |
| `backend` | str | `"matplotlib"` | Rendering backend. |
| `plink_path` | str | Auto | Path to PLINK executable. |
| `recomb_data_dir` | str | Auto | Directory with recombination maps. |
| `genomewide_threshold` | float | `5e-8` | P-value for significance line. |
| `auto_genes` | bool | `False` | If `True`, fetch the gene track with exon structure when `genes_df` is not supplied. |

### plot() Method

Create a single regional association plot. The region is passed directly;
every other option is a field of one of four frozen config models, so each is
declared once and a model built in a notebook can serve many calls.

```python
from pylocuszoom import ColumnConfig, DisplayConfig, LDConfig, PanelInputs

fig = plotter.plot(
    gwas_df,                    # Required: GWAS results
    chrom=1,                    # Required: chromosome
    start=1000000,              # Required: start position
    end=2000000,                # Required: end position
    columns=ColumnConfig(),     # Column names in gwas_df
    display=DisplayConfig(),    # Labels, overlay, gene fetching, figure size
    ld=LDConfig(),              # Lead SNP and LD source
    panels=PanelInputs(),       # Frames for the optional panels
    significance_threshold=5e-8,  # Omit to inherit the plotter; None for no line
)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gwas_df` | DataFrame | Required | GWAS results with position and p-value columns. |
| `chrom` | int or str | Required | Chromosome. |
| `start` | int | Required | Region start position (bp, `>= 1`). |
| `end` | int | Required | Region end position (bp, `> start`). |
| `columns` | `ColumnConfig` | `ColumnConfig()` | Column names, below. |
| `display` | `DisplayConfig` | `DisplayConfig()` | Display options, below. |
| `ld` | `LDConfig` | `LDConfig()` | LD options, below. |
| `panels` | `PanelInputs` | `PanelInputs()` | Optional-panel frames, below. |
| `significance_threshold` | float or None | plotter's `genomewide_threshold` | P-value for the significance line; `None` draws none. |
| `liftover` | `LiftoverConfig` | `LiftoverConfig()` | Lift source-build sumstats to the plotter's build, below. |

> **Note:** all arguments after `gwas_df` are keyword-only.

#### ColumnConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chrom_col` | str or None | `"chr"` | Chromosome column name in gwas_df. A frame without it raises; `None` selects the region by position only, for a frame already scoped to the region's chromosome. |
| `pos_col` | str | `"pos"` | Position column name in gwas_df. |
| `p_col` | str | `"p_value"` | P-value column name in gwas_df. |
| `rs_col` | str | `"rs"` | SNP ID column name in gwas_df. |

#### GenomeWideConfig

The column contract for the genome-wide families (`plot_manhattan`, `plot_qq`,
their stacked and side-by-side variants, and `plot_miami`).

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chrom_col` | str | `"chr"` | Chromosome column name. |
| `pos_col` | str | `"pos"` | Position column name. |
| `p_col` | str | `"p_value"` | P-value column name. |
| `custom_chrom_order` | list[str] | None | Chromosome order along the axis, overriding the plotter species. The canine order places PLINK's numeric sex codes beside their letters (X, 39, XY, 41, Y, 40, MT, 42); the feline order runs A1 to F2, then X, Y, MT. |

#### GenomeWideStyle

Styling for the genome-wide families, passed as `style=`. `None` keeps the
value each method draws with today, which differs by method: chromosome ticks
are 8 pt on a genomic axis and 10 pt on a category axis, categorical points
are larger, and stacked figures use smaller panel titles and axis labels.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `palette` | sequence of colours | None | Colours cycled over the chromosomes in display order, or over the categories of a categorical Manhattan. Any matplotlib colour spec, stored as hex. For a matplotlib colormap, pass its `colors`. `None` keeps the default glasbey palette. |
| `point_size` | float | None | Marker area in matplotlib `s` units for Manhattan and QQ points. |
| `point_alpha` | float | None | Marker opacity in (0, 1]. `None` draws opaque points. |
| `title_fontsize` | int | 14 | The figure title that `title=` sets on Manhattan-QQ, stacked and Miami figures. |
| `panel_title_fontsize` | int | None | Each panel's title: the QQ λ title, and the `title` of a single `plot_manhattan` or `plot_qq`. |
| `axis_label_fontsize` | int | None | X and y axis labels. |
| `tick_label_fontsize` | int | None | Tick labels on both axes. |
| `tick_step` | int | 1 | Label every n-th chromosome or category that carries data, starting with the first. At least 1. |
| `tick_rotation` | int | None | Chromosome or category tick label rotation in degrees. |
| `chrom_gap` | int | 1,000,000 | Base pairs between one chromosome's last position and the next chromosome's first. |
| `line_style` | str | `"--"` | Linestyle of the significance and suggestive lines and the QQ diagonal: `"-"`, `"--"`, `":"` or `"-."`. |
| `line_width` | float | 1.0 | Width of the same lines. |
| `title_fontweight` | str | `"bold"` | Weight of the figure title and panel titles: `"bold"` or `"normal"`. |
| `point_edge_width` | float | None | Outline width of Manhattan and QQ points; `0` draws no outline. `None` keeps 0.1 on Manhattan points and 0.02 on QQ points. |
| `y_headroom` | float | 0.1 | Space left above the highest Manhattan point or threshold line, as a fraction of its height. |
| `manhattan_qq_width_ratio` | float | 2.5 | Width of the Manhattan panel relative to the QQ panel in a Manhattan-QQ figure. |

All three backends apply every field. Corner panel labels and Miami SNP
annotations keep their own sizes.

#### DisplayConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `snp_labels` | bool | True | Whether to label top SNPs (matplotlib only). |
| `label_top_n` | int | None | Number of top SNPs to label per panel. `None` takes the method default: 5 on `plot()`, 3 on `plot_stacked()`. |
| `show_recombination` | bool | True | Whether to show recombination overlay. |
| `auto_genes` | bool | None | Fetch the gene track with exon structure when `genes_df` is not supplied; `None` inherits the constructor setting. |
| `figsize` | tuple | `(12, 8)` | Figure size (width, height) in inches. `plot_stacked()` uses the height as a floor and grows with the panel count. |

#### LDConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `lead_pos` | int | None | Lead SNP position to highlight, inside the region (`start` to `end` inclusive). Auto-detected as the strongest in-region p-value when omitted. |
| `ld_reference_file` | str | None | PLINK fileset (without extension) for LD calculation. Requires a lead. |
| `ld_col` | str | None | Column name if LD is pre-computed in gwas_df. Mutually exclusive with `ld_reference_file`. |

#### LiftoverConfig

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `lifter` | `CoordinateLifter` | None | Object with pyliftover's `convert_coordinate`, such as `pyliftover.LiftOver`. |
| `chain_path` | str or path | None | UCSC chain from the sumstats build to the plotter's build. Mutually exclusive with `lifter`; loaded chains are cached per path. |
| `lift_recombination` | bool | False | Lift the recombination maps the plotter loads through the same chain instead of the registered one. The maps must be in the chain's source build. |

#### PanelInputs

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `genes_df` | DataFrame | None | Gene annotations for track. |
| `exons_df` | DataFrame | None | Exon annotations for gene structure; needs `chr`, `start`, `end` and `gene_name`. |
| `recomb_df` | DataFrame | None | Custom recombination rate data. |
| `eqtl` | `EqtlInput` | None | The eQTL panel: `data`, `gene` (exact match on a `gene` column), `threshold` in (0, 1] (default `1e-5`), and `chrom_col` (default `"chr"`, `None` for position only). |
| `finemapping` | `FinemappingInput` | None | The fine-mapping panel: `data` with `pos` and `pip`, `cs_col` (default `"cs"`, which may be absent; another name must exist; `None` for no credible sets), and `chrom_col`. |
| `ld_heatmap` | `LDHeatmapInput` | None | The LD heatmap panel: `matrix`, `snp_ids` (required; ids must match the GWAS SNP id column, and at least one must fall inside the region, or the call raises), `height` (> 0, default `0.25`) and `metric` (`"r2"` or `"dprime"`). |

Every frame, including the ones inside the three panel models, also accepts a
PySpark DataFrame.

All config models are frozen: build a new instance rather than changing one.
`plot()` and `plot_stacked()` combine them into one internal `PlotConfig` or
`StackedPlotConfig`, which holds the rules across models (a PLINK fileset needs
a lead, every lead lies inside the region), so an invalid combination raises
`ValidationError` before anything is drawn.

### plot_stacked() Method

Create stacked plots comparing multiple GWAS. It takes the same four config
models as `plot()`, plus one list per panel-specific value.

```python
fig = plotter.plot_stacked(
    gwas_dfs,                   # Required: list of GWAS DataFrames
    chrom=1,                    # Required: chromosome
    start=1000000,              # Required: start position
    end=2000000,                # Required: end position
    columns=ColumnConfig(),     # As on plot()
    display=DisplayConfig(),    # As on plot(); label_top_n defaults to 3
    ld=LDConfig(),              # lead_pos and ld_reference_file apply to every panel
    panels=PanelInputs(),       # As on plot()
    lead_positions=None,        # Per-panel lead SNP positions
    panel_labels=None,          # Labels for each panel
    ld_reference_files=None,    # Per-panel PLINK filesets
    significance_threshold=5e-8,  # As on plot()
    liftover=LiftoverConfig(),  # As on plot(), applied to every panel
)
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gwas_dfs` | list | Required | List of GWAS DataFrames. |
| `chrom` | int or str | Required | Chromosome. |
| `start` | int | Required | Region start position (bp). |
| `end` | int | Required | Region end position (bp). |
| `columns`, `display`, `ld`, `panels` | models | defaults | As on `plot()`. |
| `lead_positions` | list | None | Lead SNP positions, one per panel. Required with a broadcast `ld.ld_reference_file`. |
| `panel_labels` | list | None | Labels, one per panel. |
| `ld_reference_files` | list | None | PLINK filesets, one per panel, replacing the broadcast `ld.ld_reference_file`. |
| `significance_threshold` | float or None | plotter's `genomewide_threshold` | As on `plot()`. |
| `liftover` | `LiftoverConfig` | `LiftoverConfig()` | As on `plot()`, applied to every frame; the window spans the lifted SNPs of all panels. |

#### ColocConfig

The value `ColocPlotter.plot_coloc(gwas_df, eqtl_df, config=...)` takes. The
two thresholds and the title are per-call arguments of `plot_coloc`, not
fields.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `gwas_p_col` | str | `"p_gwas"` | GWAS p-value column. |
| `eqtl_p_col` | str | `"p_eqtl"` | eQTL p-value column. |
| `pos_col` | str | `"pos"` | Position column, shared by both frames. |
| `rs_col` | str or None | `"rs"` | SNP id column; the default may be absent. |
| `ld_col` | str or None | None | Pre-computed LD column in `gwas_df`, which colours the points. |
| `lead_snp` | str or None | None | SNP id of the lead to highlight and label. |
| `show_correlation` | bool | True | Show the Pearson correlation. |
| `color_by_effect` | bool | False | Colour by effect-direction agreement; needs both effect columns. |
| `gwas_effect_col`, `eqtl_effect_col` | str or None | None | Effect-size columns. |
| `h4_posterior` | float or None | None | COLOC H4 posterior probability to display, in [0, 1]. |
| `figsize` | tuple | `(8, 8)` | Figure size in inches. |

### Parameter Naming Conventions

`LDConfig.lead_pos` and `ld_reference_file` are singular: on `plot()` they
describe the one panel, on `plot_stacked()` they apply to every panel.
`lead_positions` and `ld_reference_files` are the per-panel lists, one entry
per GWAS frame, and override the broadcast value.

---

## File Loaders

pyLocusZoom includes convenience functions for loading common file formats directly into DataFrames ready for plotting.

### GWAS Loaders

| Function | Format | Description |
|----------|--------|-------------|
| `load_gwas()` | Auto | Auto-detects format from filename |
| `load_plink_assoc()` | PLINK 1.x and 2 | `.assoc`, `.assoc.linear`, `.assoc.logistic`, `.qassoc`, and PLINK 2 `--glm` output (`.glm.linear`, `.glm.logistic`), whose header line starts with `#CHROM` |
| `load_regenie()` | REGENIE | `.regenie` files |
| `load_bolt_lmm()` | BOLT-LMM | `.stats` files |
| `load_gemma()` | GEMMA | `.assoc.txt` files |
| `load_saige()` | SAIGE | SAIGE output files |
| `load_gwas_catalog()` | GWAS Catalog | Summary statistics format |

```python
from pylocuszoom import load_gwas, load_plink_assoc, load_regenie

# Auto-detect format from filename extension
gwas_df = load_gwas("results.assoc.linear")

# Or use a specific loader; rename after loading if you want other names
gwas_df = load_plink_assoc("results.assoc").rename(columns={"pos": "position"})

# REGENIE (handles LOG10P conversion automatically)
gwas_df = load_regenie("ukb_chr1.regenie")
```

### eQTL Loaders

| Function | Format | Description |
|----------|--------|-------------|
| `load_gtex_eqtl()` | GTEx | Significant variant-gene pairs |
| `load_eqtl_catalogue()` | eQTL Catalogue | Standardized eQTL format |
| `load_matrixeqtl()` | MatrixEQTL | R MatrixEQTL output |

```python
from pylocuszoom import load_gtex_eqtl, load_eqtl_catalogue

# Load GTEx data, optionally filter to specific gene
eqtl_df = load_gtex_eqtl(
    "GTEx_Analysis_v8.signif_variant_gene_pairs.txt.gz",
    gene="BRCA1",  # Filter to BRCA1 associations
)

# eQTL Catalogue format
eqtl_df = load_eqtl_catalogue("eqtl_results.tsv", gene="TP53")
```

**Output columns:** `pos`, `p_value`, `gene`, `effect_size`, and `chr` when supplied.
GTEx variant IDs provide both chromosome and absolute position. A file carrying
only relative `tss_distance` is rejected because it does not locate a variant.

`calculate_colocalization_overlap(gwas_df, eqtl_df)` matches canonical chromosome
and absolute position, without allele harmonization. Custom names use
`gwas_chrom_col`, `eqtl_chrom_col` and the existing position/p-value arguments.
If both inputs are already scoped to one chromosome and omit `chr`, pass
`common_chrom=1` explicitly. Any supplied chromosome values must agree with it.
The result has `chr`, `pos`, `p_value_gwas` and `p_value_eqtl` columns. This helper
finds significant coordinate matches; it does not perform statistical colocalization.

### Fine-mapping Loaders

| Function | Format | Description |
|----------|--------|-------------|
| `load_susie()` | SuSiE | susieR output (TSV) |
| `load_finemap()` | FINEMAP | `.snp` output file |
| `load_caviar()` | CAVIAR | `.set` output file |
| `load_polyfun()` | PolyFun | PolyFun+SuSiE output |

```python
from pylocuszoom import FinemappingInput, PanelInputs, load_finemap, load_susie

# SuSiE results (handles credible set standardization)
susie_df = load_susie("susie_results.tsv")
# Output: pos, pip, cs (credible set, 0 = not in CS), rs; no chr column

# FINEMAP results (preserves PIPs and any supplied credible-set membership)
finemap_df = load_finemap("finemap_output.snp")

# SuSiE output has no chromosome column, so select its rows by position only
fig = plotter.plot_stacked(
    [gwas_df],
    chrom=1,
    start=1e6,
    end=2e6,
    panels=PanelInputs(
        finemapping=FinemappingInput(data=susie_df, chrom_col=None)
    ),
)
```

**Output columns:** `pip`, plus `pos` and `cs` where the source supplies them.
The fine-mapping panel selects rows by chromosome, and `load_susie` emits no
`chr` column: add one, or pass `FinemappingInput(..., chrom_col=None)` for a
frame that holds only the plotted chromosome.
FINEMAP and CAVIAR loaders no longer infer credible sets from cumulative PIPs.
Supply membership from the inference method that produced your results, or plot
PIPs without set assignments. CAVIAR requires a SNP annotation merge to add
absolute positions before plotting. `cs_col` chooses the output name for supplied
membership; it does not request inference.

### Gene Annotation Loaders

| Function | Format | Description |
|----------|--------|-------------|
| `load_gtf()` | GTF/GFF3 | Standard gene annotation format |
| `load_bed()` | BED | BED4+ format |
| `load_ensembl_genes()` | Ensembl | BioMart gene export |

```python
from pylocuszoom import PanelInputs, load_bed, load_gtf

# Load genes and exons from GTF
genes_df = load_gtf("gencode.v40.annotation.gtf.gz", feature_type="gene")
exons_df = load_gtf("gencode.v40.annotation.gtf.gz", feature_type="exon")

# Load from BED file
genes_df = load_bed("genes.bed")

# Use in plot
fig = plotter.plot(
    gwas_df, chrom=1, start=1e6, end=2e6,
    panels=PanelInputs(genes_df=genes_df, exons_df=exons_df),
)
```

**Output columns:** `chr`, `start`, `end`, `gene_name`, `strand` (optional)

---

## Data Formats

### GWAS Results DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str or int | Yes | Chromosome. Regional plots select the region's rows by it; pass `ColumnConfig(chrom_col=None)` for a frame that holds only the plotted chromosome. |
| `pos` | int | Yes | Genomic position (bp, 1-based). |
| `p_value` | float | Yes | P-value (0 < p ≤ 1). |
| `rs` | str | For LD/labels | SNP identifier. |

These are the canonical column names: every `load_*` function emits them and
every plotter defaults to them, so a loaded frame plots without renaming.
Other names, including the pre-4.0 `ps` and `p_wald`, are named through
`ColumnConfig` and `GenomeWideConfig`.

Regional plots select chromosome and inclusive position bounds before choosing a
lead, scaling axes, labeling points or calculating LD. A frame without the
`chrom_col` column raises; pass `ColumnConfig(chrom_col=None)` for a frame that
already holds only the requested chromosome. In stacks, shared `LDConfig`
values apply to every panel unless a per-panel list overrides them, and every
panel computing LD from a reference fileset needs a lead position. A lead
position shared by multiple variants selects the strongest p-value at that
position, with input order breaking ties. That selected row also defines label
eligibility; nearby non-lead variants are excluded before ranking labels.
Requesting reference LD replaces an existing `R2` column in the prepared plot
data. The caller's frame is unchanged. Regional heatmaps sort SNPs and both
matrix axes together, and require distinct retained genomic positions.

Every frame of a genome-wide stack is read through the same `GenomeWideConfig`
column names, and QQ compositions and Miami hover read those same columns. A
frame in other names, such as the pre-4.0 `ps` and `p_wald`, raises until you
name them. Unselected metadata never replaces a configured role. Requested colocalization LD columns
must exist in their declared source frame. Effect columns are required in their
declared sources only when `color_by_effect=True`.

Categorical Manhattan plots render missing categories as `Uncategorised`.
An explicit category order sets priority; other observed categories append in
alphabetical order so retained observations remain visible.

```python
gwas_df = pd.DataFrame({
    "chr": [1, 1, 1],
    "pos": [1000000, 1000500, 1001000],
    "p_value": [1e-8, 1e-6, 0.05],
    "rs": ["rs123", "rs456", "rs789"],
})
```

### Genes DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str/int | Yes | Chromosome (accepts "1", "chr1", or 1). |
| `start` | int | Yes | Gene start position (bp). |
| `end` | int | Yes | Gene end position (bp). |
| `gene_name` | str | Yes | Gene symbol for display. |
| `strand` | str | No | "+" or "-" for directional arrows. |
| `assembly` | str | No | Assembly the coordinates are in. Set on frames fetched from Ensembl; ignored when plotting. |

### Exons DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str/int | Yes | Chromosome. |
| `start` | int | Yes | Exon start position. |
| `end` | int | Yes | Exon end position. |
| `gene_name` | str | Yes | Parent gene (must match genes_df). |
| `assembly` | str | No | Assembly the coordinates are in. Set on frames fetched from Ensembl; ignored when plotting. |

### Recombination DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `pos` | int | Yes | Position (bp). |
| `rate` | float | Yes | Recombination rate (cM/Mb). |

### Recombination Map Files

A `recomb_data_dir` holds one tab-separated file per chromosome, named
`chr{N}_recomb.tsv` (`chr1_recomb.tsv`, `chrX_recomb.tsv`), with a header row:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome number (without "chr" prefix) |
| `pos` | Position in base pairs |
| `rate` | Recombination rate (cM/Mb) |
| `cM` | Cumulative genetic distance (optional, not used for plotting) |

```text
chr     pos     rate    cM
1       10000   0.5     0.005
1       20000   1.2     0.017
1       30000   0.8     0.025
```

### Fine-mapping DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str or int | Yes, unless `FinemappingInput(chrom_col=None)` | Chromosome. |
| `pos` | int | Yes | Variant position. |
| `pip` | float | Yes | Posterior inclusion probability (0-1). |
| `cs` | int | No | Credible set assignment (0 = not in CS). |

### eQTL DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `chr` | str or int | Yes, unless `EqtlInput(chrom_col=None)` | Chromosome. |
| `pos` | int | Yes | Variant position. |
| `p_value` | float | Yes | Association p-value. |
| `gene` | str | With `EqtlInput(gene=...)` | Target gene symbol the panel filters on. |
| `effect_size` | float | No | Effect size for color coding. |

### PheWAS DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `phenotype` | str | Yes | Phenotype name. |
| `p_value` | float | Yes | Association p-value. |
| `category` | str | No | Phenotype category for grouping/coloring. |
| `effect_size` | float | No | Effect size for direction markers. |

```python
phewas_df = pd.DataFrame({
    "phenotype": ["Height", "BMI", "T2D"],
    "p_value": [1e-15, 0.05, 1e-8],
    "category": ["Anthropometric", "Anthropometric", "Metabolic"],
})
```

### Forest Plot DataFrame

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `study` | str | Yes | Study or phenotype name. |
| `effect` | float | Yes | Effect size (beta, OR, HR). |
| `ci_lower` | float | Yes | Lower confidence interval bound. |
| `ci_upper` | float | Yes | Upper confidence interval bound. |
| `weight` | float | No | Study weight (affects marker size). |

```python
forest_df = pd.DataFrame({
    "study": ["Study A", "Study B", "Meta-analysis"],
    "effect": [0.45, 0.52, 0.46],
    "ci_lower": [0.30, 0.35, 0.40],
    "ci_upper": [0.60, 0.69, 0.52],
    "weight": [25, 35, 100],
})
```

---

## Species Support

| Species | Recombination | LD Flag | Builds |
|---------|--------------|---------|--------|
| Canine | Built-in (auto-download) | `--dog` | CanFam3.1, CanFam4 |
| Feline | User-provided | `--chr-set 18` | FelCat9 |
| Custom | User-provided | User config | Any |

### Canine

```python
# CanFam3.1 (default)
plotter = LocusZoomPlotter(species="canine")

# CanFam4 with automatic liftover
plotter = LocusZoomPlotter(species="canine", genome_build="canfam4")
```

Recombination maps from [Campbell et al. 2016](https://github.com/cflerin/dog_recombination)
are automatically downloaded on first use (~50MB), into
`recombination_maps` under the platform cache. The CanFam3.1 to CanFam4
liftover chain downloads into a `liftover` directory beside it, so replacing a
map set never touches the chain; [CONFIGURATION.md](CONFIGURATION.md#cache-location)
lists every cache folder. Managed maps requested in another assembly
without a registered conversion are skipped with a build-unavailable warning.

An explicit `recomb_data_dir` is caller-owned, read-only and already in the
requested build. This works for every species, including `species=None`, and
requires only the chromosomes used by the plot. No automatic liftover runs on
caller maps.

### Feline

Canine is the only species with built-in maps. Any other species plots its
recombination overlay from data you supply, either per plot with `recomb_df` or
from a directory of `chr{N}_recomb.tsv` files with `recomb_data_dir`.

```python
plotter = LocusZoomPlotter(species="feline")

# Provide recombination data per-plot
fig = plotter.plot(
    gwas_df, chrom=1, start=1e6, end=2e6,
    panels=PanelInputs(recomb_df=my_feline_recomb_df),
)
```

### Custom Species

```python
plotter = LocusZoomPlotter(
    species=None,
    recomb_data_dir="/path/to/recomb_maps/",
)
```

Recombination maps must be named `chr{N}_recomb.tsv`.

### Automatic Gene Annotations from Ensembl

Instead of providing your own `genes_df`, enable automatic fetching from Ensembl:

```python
plotter = LocusZoomPlotter(species="human", auto_genes=True)
fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
# Gene track populated automatically from Ensembl
```

Genes come back with their exons, so the automatic track draws intron and exon
structure rather than a plain rectangle per gene. Passing your own `exons_df`
overrides the fetched one.

**Supported Species:**

| Alias | Ensembl Name |
|-------|--------------|
| human | homo_sapiens |
| mouse | mus_musculus |
| rat | rattus_norvegicus |
| canine, dog | canis_lupus_familiaris |
| feline, cat | felis_catus |

Any valid Ensembl species name also works for annotation (e.g., `sus_scrofa`
for pig). Registered Ensembl names are aliases of their species records, so
`canis_lupus_familiaris` receives the same PLINK flags as `canine`.

LD calculation requires known chromosome-set flags. Unknown PLINK support raises
an error; supply a `Species` record with explicit `plink_flags` when adding a
species. An empty tuple explicitly selects PLINK's human defaults, while `None`
means support is unknown. Relative reference, working-directory and executable
paths resolve against the caller's directory before PLINK starts. Bare executable
names search PATH.

**Region Limit:** Maximum 5Mb per request (Ensembl API limitation). For larger regions, provide `genes_df` directly.

**Genome Build:** Genes are fetched in the build you set, from whichever source can serve it. CanFam3.1, CanFam4 and FelCat9 come from UCSC's `ncbiRefSeq` track; every other build comes from Ensembl.

| `genome_build` | Source | Assembly returned |
|----------------|--------|-------------------|
| `canfam3.1` (canine default) | UCSC `ncbiRefSeq` | CanFam3.1 |
| `canfam4`, `UU_Cfam_GSD_1.0` | UCSC `ncbiRefSeq` | UU_Cfam_GSD_1.0 |
| `felCat9` (feline default) | UCSC `ncbiRefSeq` | Felis_catus_9.0 |
| anything else | Ensembl REST | Ensembl's current assembly for the species |

Ensembl serves exactly one reference assembly per species and answers a request naming any other with that same assembly and an HTTP 200, so it cannot supply the three retired builds and will not say so. Its dog is `ROS_Cfam_1.0` and its cat is `F.catus_Fca126_mat1.0`. Release 116 was the last on the legacy REST platform and the archive REST hosts redirect to a help page, so those builds have no Ensembl source at any URL.

Both sources return the same columns, including an `assembly` column naming the assembly each row is in. On the Ensembl path, a `genome_build` that disagrees with what Ensembl served warns with a `UserWarning` naming both assemblies, since the coordinates will be off by hundreds of kilobases with nothing in the figure to show it.

UCSC's `ncbiRefSeq` is a transcript-level track, so transcripts sharing a symbol are collapsed into one gene row spanning the widest of them, and the default `biotype="protein_coding"` filter keeps genes with at least one `NM_`/`XM_` transcript. UCSC imposes no 5Mb region limit, so that cap applies only to the Ensembl path.

`genome_build` also selects the recombination map, where CanFam3.1 and CanFam4 are both supported.

**Error Handling:** A source failure raises `EnsemblAPIError` or `UCSCAPIError`, both catchable as `ReferenceAPIError`. Under `auto_genes=True` the plotter catches it, warns, and draws the plot without the gene track.

**Cache Location:** each source caches under its own folder of the platform
cache, `ensembl/{ensembl_species}/` or `ucsc/{ucsc_genome}/`, so a region
fetched from Ensembl and the same region fetched from UCSC never collide; see
[CONFIGURATION.md](CONFIGURATION.md#cache-location) for the base directory on
each platform. A CanFam3.1 or FelCat9 plot caches under `ucsc/canFam3/` or
`ucsc/felCat9/`, not under `ensembl/`. Each entry atomically publishes one ZIP containing both gene and exon CSVs.
Older separate CSV pairs become cache misses and are fetched again; clearing the
cache removes both formats. The returned count is files removed, one per new entry.

```python
# Clear cache when needed
from pylocuszoom import clear_gene_cache, get_ensembl_species_name

clear_gene_cache("ensembl")  # Every species cached from Ensembl
clear_gene_cache("ensembl", cache_species=get_ensembl_species_name("human"))
clear_gene_cache("ucsc", cache_species="canFam3")
```

The species subdirectory is named the way the source names it, so an Ensembl
one is `homo_sapiens` rather than `human`; `get_ensembl_species_name` resolves
an alias to it.

**Fetching genes without plotting:** `get_genes_for_build` returns the genes
and the exons for a region, going through the same cache the plotter uses.
`source_for` turns a species and a build into the source that can serve them,
and is the one place either is interpreted.

```python
from pylocuszoom import get_genes_for_build, source_for

genes_df, exons_df = get_genes_for_build(
    source_for("canine", "canfam3.1"), chrom=1, start=1_000_000, end=2_000_000
)
```

The result is a named tuple, so `annotations.genes` and `annotations.exons`
also work. Both frames are cached together, so a second call for the same
region makes no request at all.

**Note:** Recombination rates are NOT available from Ensembl for most species. Continue to provide recombination maps separately.

---

## Recipes & Examples

### Plot Without LD (No PLINK)

```python
fig = plotter.plot(
    gwas_df,
    chrom=1, start=1e6, end=2e6,
    # No LDConfig: neither a lead nor an LD source
)
```

All SNPs will be gray.

### Pre-computed LD

```python
from pylocuszoom import LDConfig

# LD already in the DataFrame: each SNP's r² with the lead
gwas_df = pd.DataFrame({
    "chr": [1] * 5,
    "pos": [1000000, 1000500, 1001000, 1001500, 1002000],
    "p_value": [0.05, 1e-4, 1e-8, 1e-6, 0.01],
    "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
    "r2": [0.2, 0.6, 1.0, 0.8, 0.1],
})

fig = plotter.plot(
    gwas_df,
    chrom=1, start=999_000, end=1_003_000,
    ld=LDConfig(lead_pos=1001000, ld_col="r2"),  # Use pre-computed LD
)
```

### Summary Statistics on Another Build

Summary statistics stay on the build they were computed on. To draw them over
another build's genes, give the plotter the target build and `plot()` a chain
from the sumstats build to it. `gwas_df`, `start`, `end` and `ld.lead_pos` are
then source-build coordinates; gene, eQTL and fine-mapping frames are
target-build. The window keeps the requested margins around the outermost
lifted SNPs. `plot_stacked()` takes the same `liftover` and lifts every frame;
its window spans the lifted SNPs of all panels. The config is validated in
source-build coordinates before anything is lifted, and a lead that does not
lift is an error when `ld_reference_file` needs it.

```python
from pylocuszoom import LDConfig, LiftoverConfig, LocusZoomPlotter

plotter = LocusZoomPlotter(species="canine", genome_build="canfam4")
fig = plotter.plot(
    canfam3_gwas_df,
    chrom=12, start=33_000_000, end=34_000_000,
    ld=LDConfig(lead_pos=33_500_000, ld_col="R2"),
    liftover=LiftoverConfig(chain_path="canFam3ToCanFam4.over.chain.gz"),
)
```

pyliftover positions are 0-based; the lift converts to and from the 1-based
positions in `gwas_df`. A SNP is kept only when it lifts to exactly one locus on
the same chromosome, and the lead highlight follows it. PLINK's numeric X codes
(canine 39 and 41, feline 19 and 21) and `XY` are queried as `chrX`. To see
what was dropped, call `liftover_region` yourself:

```python
from pyliftover import LiftOver
from pylocuszoom import liftover_region

lift = liftover_region(
    region_df, chrom=12, lifter=LiftOver("canFam3ToCanFam4.over.chain.gz"),
    lead_pos=33_500_000, build=plotter.genome_build,
)
lift.n_lifted, lift.n_unmapped, lift.n_multimapped, lift.n_cross_chrom
```

### Save in Multiple Formats

```python
# High-resolution PNG
fig.savefig("plot.png", dpi=300, bbox_inches="tight")

# Vector PDF for publication
fig.savefig("plot.pdf", bbox_inches="tight")

# SVG for web
fig.savefig("plot.svg", bbox_inches="tight")
```

### Logging

Logging is off by default. Turn it on, at a level, with `enable_logging`:

```python
from pylocuszoom import enable_logging

enable_logging("DEBUG")
```

`enable_logging` adds its own stderr sink, and pyLocusZoom's records also reach
every sink without a filter. In a script that has not touched loguru, its
default stderr sink is still installed, so each line prints twice: once in
loguru's default format and once as `LEVEL | pylocuszoom | message`. Call
`logger.remove(0)` from loguru first to keep only the second, or skip
`enable_logging` and call `logger.enable("pylocuszoom")` to route the records
through your own sinks.

### Custom Significance Threshold

`genomewide_threshold` sets where that plotter draws its significance line.
`LocusZoomPlotter`, `ManhattanPlotter`, `MiamiPlotter`, `StatsPlotter`, and
`ColocPlotter` all take it. `ColocPlotter` also takes `eqtl_threshold` for its
horizontal line. It applies to the plot families that have a significance line; QQ plots
(`plot_qq`) and forest plots (`plot_forest`) do not draw one and ignore it.

```python
# Suggestive threshold
plotter = LocusZoomPlotter(
    species="canine",
    genomewide_threshold=1e-5,
)
```

Every plot method with a significance line also takes a per-call threshold
argument. Omit it to inherit the plotter's `genomewide_threshold`, pass a
p-value to override it for one call, or pass `None` to draw no line. The
argument is not spelled the same way everywhere, and two methods do not take
one at all:

| Method | Per-call argument |
|--------|-------------------|
| `LocusZoomPlotter.plot` | `significance_threshold` |
| `LocusZoomPlotter.plot_stacked` | `significance_threshold` |
| `ManhattanPlotter.plot_manhattan` | `significance_threshold` |
| `ManhattanPlotter.plot_manhattan_stacked` | `significance_threshold` |
| `ManhattanPlotter.plot_manhattan_qq` | `significance_threshold` |
| `ManhattanPlotter.plot_manhattan_qq_stacked` | `significance_threshold` |
| `StatsPlotter.plot_phewas` | `significance_threshold` |
| `MiamiPlotter.plot_miami` | `top_threshold` and `bottom_threshold` |
| `ColocPlotter.plot_coloc` | `gwas_threshold` and `eqtl_threshold` |
| `ManhattanPlotter.plot_qq` | none |
| `StatsPlotter.plot_forest` | none |

```python
plotter = ManhattanPlotter(genomewide_threshold=1e-5)

plotter.plot_manhattan(df)                                # line at 1e-5
plotter.plot_manhattan(df, significance_threshold=5e-8)   # line at 5e-8
plotter.plot_manhattan(df, significance_threshold=None)   # no line

# Miami names one threshold per panel, not significance_threshold
miami = MiamiPlotter(genomewide_threshold=1e-5)
miami.plot_miami(top_df, bottom_df, top_threshold=5e-8, bottom_threshold=None)
```

> **Changed:** before this release, `ManhattanPlotter`, `MiamiPlotter`, and
> `StatsPlotter` accepted `genomewide_threshold` and then ignored it, always
> drawing at 5e-8. If you passed it and worked around the old behaviour by also
> passing the per-call argument, that still works and still wins.
> `ColocPlotter` had no constructor threshold at all and its two per-call
> arguments were plain floats, so there was no way to ask it for no line.

### Large Datasets with PySpark

```python
from pylocuszoom import LocusZoomPlotter, to_pandas

# Automatic conversion
fig = plotter.plot(spark_df, chrom=1, start=1e6, end=2e6)

# Manual conversion with sampling
pandas_df = to_pandas(spark_df, sample_size=100000)
fig = plotter.plot(pandas_df, chrom=1, start=1e6, end=2e6)
```

---

## Troubleshooting

### PLINK Not Found

```text
PlinkError: PLINK not found. Install PLINK 1.9 or specify plink_path.
```

Install PLINK 1.9 and add to PATH, or specify:

```python
plotter = LocusZoomPlotter(plink_path="/path/to/plink")
```

### Missing Recombination Maps

The first canine plot triggers an automatic download. If it fails:

```python
from pylocuszoom import download_canine_recombination_maps
download_canine_recombination_maps()
```

### Figures in Jupyter

pyLocusZoom changes no notebook display setting. A matplotlib figure is a
pyplot figure, so the inline backend shows it when the cell ends; ending the
cell with `fig` as well shows it a second time. A plotly figure shows as the
cell's last expression or through `fig.show()`. A bokeh figure shows through
`bokeh.io.show(fig)` after one `bokeh.io.output_notebook()` call.

### LD Calculation Fails

When LD cannot colour a panel (PLINK finds no pairs for the lead, or the
reference panel names a variant twice), the plot is drawn without LD colouring
and one `UserWarning` names the panel and the reason. A GWAS frame without its
SNP-id column raises `ValidationError` before any PLINK call, and a PLINK
failure (not found, non-zero exit, timeout) raises `PlinkError`. Ensure:

1. GWAS DataFrame has `rs` column (or specify `rs_col`)
2. SNP IDs match those in PLINK fileset
3. Lead SNP exists in both datasets

---

## API Stability

`pylocuszoom.__all__` exports two tiers. Both work and both are tested; they
differ in what a minor release is allowed to do to them. `__all__` in
`pylocuszoom/__init__.py` is written in the same two blocks.

### Core

The plot families, the values you hand them, the loaders that produce those
values, and the errors they raise. These follow semantic versioning: a break
here means a major release, with a CHANGELOG entry and a migration note.

| Group | Names |
|-------|-------|
| Plotters | `LocusZoomPlotter`, `ManhattanPlotter`, `MiamiPlotter`, `StatsPlotter`, `LDHeatmapPlotter`, `ColocPlotter` |
| Plot configuration | `ColumnConfig`, `DisplayConfig`, `GenomeWideConfig`, `GenomeWideStyle`, `LDConfig`, `LiftoverConfig`, `PanelInputs`, `EqtlInput`, `FinemappingInput`, `LDHeatmapInput`, `ColocConfig` |
| Column vocabulary | `Canonical` |
| GWAS loaders | `load_gwas`, `load_plink_assoc`, `load_regenie`, `load_bolt_lmm`, `load_gemma`, `load_saige`, `load_gwas_catalog` |
| eQTL loaders | `load_gtex_eqtl`, `load_eqtl_catalogue`, `load_matrixeqtl` |
| Fine-mapping loaders | `load_susie`, `load_finemap`, `load_caviar`, `load_polyfun` |
| Gene annotation loaders | `load_gtf`, `load_bed`, `load_ensembl_genes` |
| Species | `Species`, `resolve_species` |
| Logging | `enable_logging`, `disable_logging` |
| Exceptions | `PyLocusZoomError`, `ValidationError`, `DataDownloadError`, `EmptyLDOutputError`, `EnsemblAPIError`, `OptionalDependencyMissing`, `ReferenceAPIError`, `UCSCAPIError`, `PlinkError`, `PheWASValidationError`, `ForestValidationError`, `EQTLValidationError`, `FinemappingValidationError`, `LoaderValidationError`, `LDUnavailableError`, `RecombinationMapNotFound` |
| Metadata | `__version__` |

### Toolbox

Internals that are useful on their own and exported because somebody wanted
them. None is deprecated and none is dead, but any of them may change shape or
move between minor releases. Building on one is fine; pin the version if you
do, and open an issue so the name can be promoted to core.

| Group | Names |
|-------|-------|
| Backends | `BackendType`, `get_backend` |
| Colours | `get_ld_color`, `get_ld_bin`, `get_ld_color_palette`, `get_phewas_category_color`, `get_phewas_category_palette`, `LEAD_SNP_COLOR`, `PHEWAS_CATEGORY_COLORS` |
| LD | `calculate_ld`, `calculate_pairwise_ld` |
| SNP labels | `add_snp_labels`, `adjust_snp_labels` |
| Gene track | `get_nearest_gene` |
| Gene reference routing | `get_genes_for_build`, `source_for`, `clear_gene_cache`, `get_ensembl_species_name` |
| Recombination maps | `download_canine_recombination_maps`, `ensure_recomb_maps`, `get_recombination_rate_for_region`, `load_recombination_map` |
| Liftover | `CoordinateLifter`, `liftover_region`, `RegionLiftResult` |
| eQTL helpers | `filter_eqtl_by_gene`, `filter_eqtl_by_region`, `prepare_eqtl_for_plotting`, `get_eqtl_genes`, `calculate_colocalization_overlap` |
| Fine-mapping helpers | `filter_finemapping_by_region`, `filter_by_credible_set`, `get_credible_sets`, `get_top_pip_variants`, `prepare_finemapping_for_plotting` |
| DataFrame helpers | `to_pandas` |

---

## See Also

- [Code Map](CODEMAP.md) - Architecture and source code navigation
- [Example Notebook](../examples/getting_started.ipynb) - Interactive tutorial
- [GitHub Issues](https://github.com/michael-denyer/pylocuszoom/issues) - Bug reports
- [CHANGELOG](../CHANGELOG.md) - Version history

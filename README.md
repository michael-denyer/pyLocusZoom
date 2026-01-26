# snp-scope-plot

<img src="logo.svg" alt="snp-scope-plot logo" width="120" align="right">

Regional association plots for GWAS results with LD coloring, gene tracks, and recombination rate overlays.

Inspired by [LocusZoom](http://locuszoom.org/) and [LocusZoomR](https://github.com/Kungsgansen/locuszoomr).

## Features

- **LD coloring**: SNPs colored by linkage disequilibrium (R²) with lead variant
- **Gene track**: Annotated gene/exon positions below the association plot
- **Recombination rate**: Overlay showing recombination rate across region (dog only)
- **SNP labels**: Automatic labeling of top SNPs with RS ID or nearest gene
- **Species support**: Built-in dog (CanFam3.1/CanFam4), cat (FelCat9), or custom species
- **CanFam4 support**: Automatic coordinate liftover for recombination maps

## Installation

```bash
pip install snp-scope-plot
```

Or install from source:

```bash
pip install -e packages/snp-scope-plot
```

## Quick Start

```python
from snp_scope_plot import SNPScopePlotter

# Initialize plotter (loads reference data for dog)
plotter = SNPScopePlotter(species="dog")

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
from snp_scope_plot import SNPScopePlotter

plotter = SNPScopePlotter(
    species="dog",                      # or "cat", or None for custom
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

The default genome build for dog is CanFam3.1. For CanFam4 data:

```python
plotter = SNPScopePlotter(species="dog", genome_build="canfam4")
```

Recombination maps are automatically lifted over from CanFam3.1 to CanFam4 coordinates using the UCSC liftOver chain file.

## Using with Other Species

```python
# Cat (LD and gene tracks, user provides recombination data)
plotter = SNPScopePlotter(species="cat")

# Custom species (provide all reference data)
plotter = SNPScopePlotter(
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

Dog recombination maps are downloaded from [Campbell et al. 2016](https://github.com/cflerin/dog_recombination) on first use.

To manually download:

```python
from snp_scope_plot import download_dog_recombination_maps

download_dog_recombination_maps()
```

## Logging

Logging uses [loguru](https://github.com/Delgan/loguru) and is configured via the `log_level` parameter (default: `"INFO"`):

```python
# Suppress logging
plotter = SNPScopePlotter(log_level=None)

# Enable DEBUG level for troubleshooting
plotter = SNPScopePlotter(log_level="DEBUG")
```

## Requirements

- Python >= 3.9
- matplotlib >= 3.5.0
- pandas >= 1.4.0
- numpy >= 1.21.0
- loguru >= 0.7.0
- pyliftover >= 0.4 (for CanFam4 coordinate liftover)
- [PLINK 1.9](https://www.cog-genomics.org/plink/) (for LD calculations) - must be on PATH or specify `plink_path`

## License

MIT

# pyLocusZoom Architecture

## Project Structure

```
pyLocusZoom/
├── .github/workflows/
│   ├── ci.yml                    # CI pipeline (tests, lint)
│   └── publish.yml               # PyPI publish (Trusted Publishing)
│
├── src/pylocuszoom/
│   ├── __init__.py          (120) # Public API exports
│   ├── plotter.py           (733) # LocusZoomPlotter - main entry point
│   │                              #   plot(), plot_stacked()
│   ├── backends/
│   │   ├── __init__.py       (52) # Backend registry, get_backend()
│   │   ├── base.py          (341) # PlotBackend protocol
│   │   ├── matplotlib_backend.py (288) # Static plots
│   │   ├── plotly_backend.py    (474) # Interactive with hover
│   │   └── bokeh_backend.py     (441) # Dashboard-friendly
│   │
│   ├── colors.py            (107) # LD color palettes
│   ├── gene_track.py        (311) # Gene/exon rendering
│   ├── labels.py            (118) # SNP label positioning
│   ├── ld.py                (209) # PLINK LD calculation
│   ├── eqtl.py              (218) # eQTL data handling
│   ├── recombination.py     (432) # Recomb map loading/liftover
│   ├── logging.py           (153) # Loguru configuration
│   ├── utils.py             (194) # Validation, PySpark support
│   └── reference_data/        (4) # Cached recomb maps location
│
├── tests/                   (~1.2k) # pytest suite
├── examples/
│   └── getting_started.ipynb      # Tutorial notebook
│
├── pyproject.toml                 # Build config, dependencies
├── README.md                      # Documentation
└── LICENSE.md                     # GPL-3.0-or-later
```

## Architecture Diagram

```mermaid
graph TD
    subgraph User API
        LZP[LocusZoomPlotter]
    end

    subgraph Backends
        MPL[MatplotlibBackend]
        PLY[PlotlyBackend]
        BOK[BokehBackend]
    end

    subgraph Components
        LD[ld.py<br/>calculate_ld]
        GT[gene_track.py<br/>plot_gene_track]
        RC[recombination.py<br/>add_recombination_overlay]
        LB[labels.py<br/>add_snp_labels]
        EQ[eqtl.py<br/>prepare_eqtl_for_plotting]
        CO[colors.py<br/>get_ld_color]
    end

    LZP --> MPL
    LZP --> PLY
    LZP --> BOK

    LZP --> LD
    LZP --> GT
    LZP --> RC
    LZP --> LB
    LZP --> EQ
    LZP --> CO
```

## Data Flow

```mermaid
flowchart LR
    subgraph Input
        GWAS[GWAS DataFrame]
        GENES[Genes DataFrame]
        EQTL[eQTL DataFrame]
        PLINK[PLINK files]
    end

    subgraph Processing
        VAL[Validation]
        LDC[LD Calculation]
        REC[Recombination Lookup]
        COL[Color Assignment]
    end

    subgraph Rendering
        BE[Backend]
        SC[Scatter Plot]
        GT[Gene Track]
        OV[Recomb Overlay]
        LB[Labels]
    end

    subgraph Output
        FIG[Figure]
    end

    GWAS --> VAL --> LDC
    PLINK --> LDC
    LDC --> COL --> SC
    GENES --> GT
    EQTL --> SC
    VAL --> REC --> OV
    SC --> BE
    GT --> BE
    OV --> BE
    LB --> BE
    BE --> FIG
```

## Key Entry Points

| Function | Location | Purpose |
|----------|----------|---------|
| `LocusZoomPlotter()` | `plotter.py:54` | Main constructor |
| `.plot()` | `plotter.py:200` | Single regional plot |
| `.plot_stacked()` | `plotter.py:450` | Multi-GWAS stacked plot |
| `get_backend()` | `backends/__init__.py:30` | Backend factory |

## Backend Protocol

All backends implement the `PlotBackend` protocol defined in `backends/base.py`:

```python
class PlotBackend(Protocol):
    def create_figure(self, n_panels, height_ratios, figsize) -> Any
    def create_scatter(self, ax, x, y, colors, sizes, hover_data, **kwargs) -> None
    def create_line(self, ax, x, y, **kwargs) -> None
    def add_annotation(self, ax, text, xy, **kwargs) -> None
    def finalize(self, fig, **kwargs) -> Any
    def save(self, fig, path, **kwargs) -> None
```

## Module Responsibilities

| Module | Responsibility |
|--------|----------------|
| `plotter.py` | Orchestrates plot creation, manages backends |
| `backends/` | Rendering abstraction for matplotlib/plotly/bokeh |
| `ld.py` | PLINK subprocess for LD calculation |
| `gene_track.py` | Gene/exon rectangle rendering |
| `recombination.py` | Load/cache/liftover recombination maps |
| `colors.py` | LD-to-color mapping with standard palette |
| `labels.py` | Non-overlapping SNP label placement |
| `eqtl.py` | eQTL data validation and preparation |
| `utils.py` | DataFrame validation, PySpark conversion |
| `logging.py` | Loguru configuration |

## Dependencies

### Required
- matplotlib >= 3.5.0
- pandas >= 1.4.0
- numpy >= 1.21.0
- loguru >= 0.7.0
- plotly >= 5.0.0
- bokeh >= 3.8.2
- kaleido >= 0.2.0
- pyliftover >= 0.4

### Optional
- pyspark >= 3.0.0 (for large-scale data)
- adjustText >= 0.8 (for improved label placement)

### External
- PLINK 1.9 (for LD calculations)

# pyLocusZoom Code Map

Bidirectional navigation between architecture diagram and source code.

## System Overview

```mermaid
flowchart TB
    subgraph Input["Input Layer [1]"]
        A["LocusZoomPlotter [1a]"]
        B["PlotConfig [1b]"]
        C["File Loaders [1c]"]
    end

    subgraph Validation["Validation Layer [2]"]
        D["DataFrameValidator [2a]"]
        E["Pydantic Schemas [2b]"]
    end

    subgraph Core["Core Processing [3]"]
        F["LD Calculator [3a]"]
        G["Color Assignment [3b]"]
        H["Gene Layout [3c]"]
        I["Recombination [3d]"]
    end

    subgraph Backends["Backend Layer [4]"]
        J["PlotBackend Protocol [4a]"]
        K["MatplotlibBackend [4b]"]
        L["PlotlyBackend [4c]"]
        M["BokehBackend [4d]"]
        N["HoverDataBuilder [4e]"]
    end

    subgraph Output["Output [5]"]
        O["Static PNG/PDF"]
        P["Interactive HTML"]
    end

    A --> B
    C --> D
    B --> D
    D --> E
    E --> F
    E --> G
    E --> H
    E --> I
    F --> J
    G --> J
    H --> J
    I --> J
    J --> K
    J --> L
    J --> M
    L --> N
    M --> N
    K --> O
    L --> P
    M --> P
```

---

## [1] Input Layer (Entry Points)

User-facing API for creating plots.

| ID | Component | Description | File:Line |
|----|-----------|-------------|-----------|
| 1a | LocusZoomPlotter | Main plotter class with plot methods | [plotter.py:62](../src/pylocuszoom/plotter.py#L62) |
| 1b | PlotConfig | Pydantic configuration for plot() | [config.py:122](../src/pylocuszoom/config.py#L122) |
| 1b | StackedPlotConfig | Configuration for plot_stacked() | [config.py:233](../src/pylocuszoom/config.py#L233) |
| 1c | load_gwas | Auto-detect GWAS format loader | [loaders.py:45](../src/pylocuszoom/loaders.py#L45) |
| 1c | load_susie | SuSiE fine-mapping loader | [loaders.py:520](../src/pylocuszoom/loaders.py#L520) |
| 1c | load_gtex_eqtl | GTEx eQTL loader | [loaders.py:380](../src/pylocuszoom/loaders.py#L380) |

### Config Components [1b]

| Config | Purpose | File:Line |
|--------|---------|-----------|
| RegionConfig | chrom, start, end validation | [config.py:21](../src/pylocuszoom/config.py#L21) |
| ColumnConfig | pos_col, p_col, rs_col mapping | [config.py:44](../src/pylocuszoom/config.py#L44) |
| DisplayConfig | snp_labels, figsize, show_recomb | [config.py:60](../src/pylocuszoom/config.py#L60) |
| LDConfig | lead_pos, ld_reference_file | [config.py:82](../src/pylocuszoom/config.py#L82) |

### File Loaders [1c]

| Module | Formats | File |
|--------|---------|------|
| GWAS | PLINK, REGENIE, BOLT-LMM, GEMMA, SAIGE | [loaders.py](../src/pylocuszoom/loaders.py) |
| eQTL | GTEx, eQTL Catalogue, MatrixEQTL | [loaders.py](../src/pylocuszoom/loaders.py) |
| Fine-mapping | SuSiE, FINEMAP, CAVIAR, PolyFun | [loaders.py](../src/pylocuszoom/loaders.py) |
| Genes | GTF, BED, Ensembl | [loaders.py](../src/pylocuszoom/loaders.py) |

---

## [2] Validation Layer

Input validation and schema enforcement.

| ID | Component | Description | File:Line |
|----|-----------|-------------|-----------|
| 2a | DataFrameValidator | Fluent builder for DataFrame checks | [validation.py:15](../src/pylocuszoom/validation.py#L15) |
| 2b | GWASRowModel | Pydantic model for GWAS rows | [schemas.py:20](../src/pylocuszoom/schemas.py#L20) |
| 2b | EQTLRowModel | Pydantic model for eQTL rows | [schemas.py:65](../src/pylocuszoom/schemas.py#L65) |
| 2b | FinemappingRowModel | Pydantic model for fine-mapping | [schemas.py:95](../src/pylocuszoom/schemas.py#L95) |

### DataFrameValidator Methods [2a]

```python
validator = DataFrameValidator(df, "gwas_df")
validator.require_columns(["chr", "pos", "p"])
        .require_numeric(["pos", "p"])
        .require_range("p", min_val=0, max_val=1)
        .validate()  # Raises ValidationError
```

---

## [3] Core Processing

Data transformation and enrichment.

| ID | Component | Description | File:Line |
|----|-----------|-------------|-----------|
| 3a | calculate_ld | PLINK wrapper for LD calculation | [ld.py:120](../src/pylocuszoom/ld.py#L120) |
| 3a | find_plink | Locate PLINK executable | [ld.py:25](../src/pylocuszoom/ld.py#L25) |
| 3b | get_ld_color | Map R² to hex color | [colors.py:45](../src/pylocuszoom/colors.py#L45) |
| 3b | get_credible_set_color | CS index to color | [colors.py:120](../src/pylocuszoom/colors.py#L120) |
| 3c | assign_gene_positions | Layout genes to rows | [gene_track.py:85](../src/pylocuszoom/gene_track.py#L85) |
| 3c | plot_gene_track_generic | Render gene track | [gene_track.py:180](../src/pylocuszoom/gene_track.py#L180) |
| 3d | get_recombination_rate_for_region | Load recomb data | [recombination.py:150](../src/pylocuszoom/recombination.py#L150) |
| 3d | download_canine_recombination_maps | Auto-download maps | [recombination.py:85](../src/pylocuszoom/recombination.py#L85) |

### LD Color Scheme [3b]

```python
LD_BINS = [
    (0.8, "0.8 - 1.0", "#FF0000"),  # red
    (0.6, "0.6 - 0.8", "#FFA500"),  # orange
    (0.4, "0.4 - 0.6", "#00CD00"),  # green
    (0.2, "0.2 - 0.4", "#00EEEE"),  # cyan
    (0.0, "0.0 - 0.2", "#4169E1"),  # blue
]
LEAD_SNP_COLOR = "#7D26CD"  # purple
```

---

## [4] Backend Layer

Rendering implementations.

| ID | Component | Description | File:Line |
|----|-----------|-------------|-----------|
| 4a | PlotBackend | Protocol defining backend interface | [base.py:11](../src/pylocuszoom/backends/base.py#L11) |
| 4b | MatplotlibBackend | Static publication plots | [matplotlib_backend.py:25](../src/pylocuszoom/backends/matplotlib_backend.py#L25) |
| 4c | PlotlyBackend | Interactive HTML (hover) | [plotly_backend.py:30](../src/pylocuszoom/backends/plotly_backend.py#L30) |
| 4d | BokehBackend | Dashboard integration | [bokeh_backend.py:28](../src/pylocuszoom/backends/bokeh_backend.py#L28) |
| 4e | HoverDataBuilder | Construct hover tooltips | [hover.py:25](../src/pylocuszoom/backends/hover.py#L25) |

### Backend Protocol Methods [4a]

| Method | Purpose |
|--------|---------|
| `create_figure()` | Initialize figure & axes |
| `scatter()` | Plot points with LD colors |
| `axhline()` | Significance threshold line |
| `secondary_yaxis()` | Recombination rate axis |
| `add_snp_labels()` | SNP ID labels (matplotlib) |
| `add_ld_legend()` | LD color legend |
| `finalize_layout()` | Apply margins & spacing |

### Backend Capabilities

| Backend | SNP Labels | Hover | Secondary Axis |
|---------|------------|-------|----------------|
| matplotlib | ✅ | ❌ | ✅ |
| plotly | ❌ | ✅ | ✅ |
| bokeh | ❌ | ✅ | ✅ |

---

## [5] Supporting Infrastructure

| Component | Description | File |
|-----------|-------------|------|
| Exceptions | PyLocusZoomError hierarchy | [exceptions.py](../src/pylocuszoom/exceptions.py) |
| Logging | loguru-based logger | [logging.py](../src/pylocuszoom/logging.py) |
| Ensembl API | Gene fetching from REST | [ensembl.py](../src/pylocuszoom/ensembl.py) |
| Utilities | to_pandas, normalize_chrom | [utils.py](../src/pylocuszoom/utils.py) |

### Exception Hierarchy

```
PyLocusZoomError (base)
├── ValidationError (also ValueError)
│   ├── EQTLValidationError
│   ├── FinemappingValidationError
│   └── LoaderValidationError
├── BackendError
└── DataDownloadError (also RuntimeError)
```

---

## Data Flow Diagram

```mermaid
sequenceDiagram
    participant User
    participant Plotter as LocusZoomPlotter [1a]
    participant Config as PlotConfig [1b]
    participant Valid as DataFrameValidator [2a]
    participant LD as calculate_ld [3a]
    participant Backend as PlotBackend [4a]

    User->>Plotter: plot(gwas_df, config)
    Plotter->>Config: validate configuration
    Plotter->>Valid: validate_gwas_df()
    Valid-->>Plotter: validated DataFrame

    opt LD calculation requested
        Plotter->>LD: calculate_ld(lead_pos)
        LD-->>Plotter: R² values
    end

    Plotter->>Backend: create_figure()
    Plotter->>Backend: scatter(positions, -log10p, colors)
    Plotter->>Backend: axhline(significance)

    opt Gene track
        Plotter->>Backend: plot_gene_track()
    end

    opt Recombination
        Plotter->>Backend: secondary_yaxis()
    end

    Backend-->>User: figure object
```

---

## Quick Navigation

| Area | Entry Point |
|------|-------------|
| Main API | [plotter.py](../src/pylocuszoom/plotter.py) |
| Configuration | [config.py](../src/pylocuszoom/config.py) |
| Validation | [validation.py](../src/pylocuszoom/validation.py) |
| All loaders | [loaders.py](../src/pylocuszoom/loaders.py) |
| Backend protocol | [backends/base.py](../src/pylocuszoom/backends/base.py) |
| Colors & themes | [colors.py](../src/pylocuszoom/colors.py) |
| Gene rendering | [gene_track.py](../src/pylocuszoom/gene_track.py) |
| LD calculation | [ld.py](../src/pylocuszoom/ld.py) |
| Exceptions | [exceptions.py](../src/pylocuszoom/exceptions.py) |

---

## File Statistics

| Category | Files | Lines |
|----------|-------|-------|
| Core | 15 | ~4,500 |
| Backends | 5 | ~2,800 |
| Tests | 15+ | ~3,000 |
| **Total** | **35+** | **~10,300** |

"""pyLocusZoom - Regional association plots for GWAS results.

This package provides LocusZoom-style regional association plots with:
- LD coloring based on R² with lead variant
- Gene and exon tracks
- Recombination rate overlays (canine built-in, or user-provided)
- Automatic SNP labeling
- Multiple backends: matplotlib (static), plotly (interactive), bokeh (dashboards)
- eQTL overlay support
- Fine-mapping/SuSiE visualization (PIP line with credible set coloring)
- PySpark DataFrame support for large-scale data

Example:
    >>> from pylocuszoom import LocusZoomPlotter
    >>> plotter = LocusZoomPlotter(species="canine")
    >>> fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
    >>> fig.savefig("regional_plot.png", dpi=150)

Interactive example:
    >>> plotter = LocusZoomPlotter(species="canine", backend="plotly")
    >>> fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
    >>> fig.write_html("regional_plot.html")

Stacked plots:
    >>> fig = plotter.plot_stacked(
    ...     [gwas_height, gwas_bmi],
    ...     chrom=1, start=1000000, end=2000000,
    ...     panel_labels=["Height", "BMI"],
    ... )

Species Support:
    - Canine (Canis lupus familiaris): Full features including built-in recombination maps
    - Feline (Felis catus): LD coloring and gene tracks (user provides recombination data)
    - Custom: User provides all reference data
"""

__version__ = "0.6.0"

# Main plotter class
# Backend types
from .backends import BackendType, get_backend

# Colors and LD
from .colors import (
    LEAD_SNP_COLOR,
    PHEWAS_CATEGORY_COLORS,
    get_ld_bin,
    get_ld_color,
    get_ld_color_palette,
    get_phewas_category_color,
    get_phewas_category_palette,
)

# Configuration classes (internal - use kwargs directly on plot()/plot_stacked())
# from .config import PlotConfig, StackedPlotConfig  # Internal use only
# Ensembl integration
from .ensembl import (
    clear_ensembl_cache,
    fetch_exons_from_ensembl,
    fetch_genes_from_ensembl,
    get_ensembl_species_name,
    get_genes_for_region,
)

# eQTL support
from .eqtl import (
    calculate_colocalization_overlap,
    filter_eqtl_by_gene,
    filter_eqtl_by_region,
    get_eqtl_genes,
    prepare_eqtl_for_plotting,
    validate_eqtl_df,
)

# Exception hierarchy
from .exceptions import (
    BackendError,
    DataDownloadError,
    EQTLValidationError,
    FinemappingValidationError,
    LoaderValidationError,
    PyLocusZoomError,
    ValidationError,
)

# Fine-mapping/SuSiE support
from .finemapping import (
    filter_by_credible_set,
    filter_finemapping_by_region,
    get_credible_sets,
    get_top_pip_variants,
    prepare_finemapping_for_plotting,
    validate_finemapping_df,
)

# Forest plot support
from .forest import validate_forest_df

# Gene track
from .gene_track import get_nearest_gene, plot_gene_track

# Labels
from .labels import add_snp_labels

# LD calculation
from .ld import calculate_ld

# File format loaders
from .loaders import (
    load_bed,
    load_bolt_lmm,
    load_caviar,
    load_ensembl_genes,
    load_eqtl_catalogue,
    load_finemap,
    load_gemma,
    # eQTL loaders
    load_gtex_eqtl,
    # Gene annotation loaders
    load_gtf,
    # GWAS loaders
    load_gwas,
    load_gwas_catalog,
    load_matrixeqtl,
    load_plink_assoc,
    load_polyfun,
    load_regenie,
    load_saige,
    # Fine-mapping loaders
    load_susie,
)

# Logging configuration
from .logging import disable_logging, enable_logging

# PheWAS support
from .phewas import validate_phewas_df
from .plotter import LocusZoomPlotter

# Reference data management
from .recombination import (
    add_recombination_overlay,
    download_canine_recombination_maps,
    get_recombination_rate_for_region,
    load_recombination_map,
)

# Validation utilities
from .utils import to_pandas

__all__ = [
    # Core
    "__version__",
    "LocusZoomPlotter",
    # Backends
    "BackendType",
    "get_backend",
    # Reference data
    "download_canine_recombination_maps",
    # Colors
    "get_ld_color",
    "get_ld_bin",
    "get_ld_color_palette",
    "get_phewas_category_color",
    "get_phewas_category_palette",
    "LEAD_SNP_COLOR",
    "PHEWAS_CATEGORY_COLORS",
    # Gene track
    "get_nearest_gene",
    "plot_gene_track",
    # LD
    "calculate_ld",
    # Labels
    "add_snp_labels",
    # Recombination
    "add_recombination_overlay",
    "get_recombination_rate_for_region",
    "load_recombination_map",
    # eQTL
    "validate_eqtl_df",
    "filter_eqtl_by_gene",
    "filter_eqtl_by_region",
    "prepare_eqtl_for_plotting",
    "get_eqtl_genes",
    "calculate_colocalization_overlap",
    "EQTLValidationError",
    # Fine-mapping/SuSiE
    "validate_finemapping_df",
    "filter_finemapping_by_region",
    "filter_by_credible_set",
    "get_credible_sets",
    "get_top_pip_variants",
    "prepare_finemapping_for_plotting",
    "FinemappingValidationError",
    # Logging
    "enable_logging",
    "disable_logging",
    # Exceptions
    "PyLocusZoomError",
    "ValidationError",
    "BackendError",
    "DataDownloadError",
    # Utils
    "to_pandas",
    # PheWAS
    "validate_phewas_df",
    # Forest plot
    "validate_forest_df",
    # GWAS loaders
    "load_gwas",
    "load_plink_assoc",
    "load_regenie",
    "load_bolt_lmm",
    "load_gemma",
    "load_saige",
    "load_gwas_catalog",
    # eQTL loaders
    "load_gtex_eqtl",
    "load_eqtl_catalogue",
    "load_matrixeqtl",
    # Fine-mapping loaders
    "load_susie",
    "load_finemap",
    "load_caviar",
    "load_polyfun",
    # Gene annotation loaders
    "load_gtf",
    "load_bed",
    "load_ensembl_genes",
    # Schema validation
    "LoaderValidationError",
    # Ensembl integration
    "get_genes_for_region",
    "fetch_genes_from_ensembl",
    "fetch_exons_from_ensembl",
    "get_ensembl_species_name",
    "clear_ensembl_cache",
]

"""pyLocusZoom - Regional association plots for GWAS results.

This package provides LocusZoom-style regional association plots with:
- LD coloring based on R² with lead variant
- Gene and exon tracks
- Recombination rate overlays (dog built-in, or user-provided)
- Automatic SNP labeling
- Multiple backends: matplotlib (static), plotly (interactive), bokeh (dashboards)
- eQTL overlay support
- PySpark DataFrame support for large-scale data

Example:
    >>> from pylocuszoom import LocusZoomPlotter
    >>> plotter = LocusZoomPlotter(species="dog")
    >>> fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
    >>> fig.savefig("regional_plot.png", dpi=150)

Interactive example:
    >>> plotter = LocusZoomPlotter(species="dog", backend="plotly")
    >>> fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
    >>> fig.write_html("regional_plot.html")

Stacked plots:
    >>> fig = plotter.plot_stacked(
    ...     [gwas_height, gwas_bmi],
    ...     chrom=1, start=1000000, end=2000000,
    ...     panel_labels=["Height", "BMI"],
    ... )

Species Support:
    - Dog (Canis lupus familiaris): Full features including built-in recombination maps
    - Cat (Felis catus): LD coloring and gene tracks (user provides recombination data)
    - Custom: User provides all reference data
"""

__version__ = "0.1.0"

# Main plotter class
from .plotter import LocusZoomPlotter

# Backend types
from .backends import BackendType, get_backend

# Colors and LD
from .colors import LEAD_SNP_COLOR, get_ld_bin, get_ld_color, get_ld_color_palette

# Gene track
from .gene_track import get_nearest_gene, plot_gene_track

# Labels
from .labels import add_snp_labels

# LD calculation
from .ld import calculate_ld

# Logging configuration
from .logging import disable_logging, enable_logging

# Reference data management
from .recombination import (
    add_recombination_overlay,
    download_dog_recombination_maps,
    get_recombination_rate_for_region,
    load_recombination_map,
)

# eQTL support
from .eqtl import (
    EQTLValidationError,
    calculate_colocalization_overlap,
    filter_eqtl_by_gene,
    filter_eqtl_by_region,
    get_eqtl_genes,
    prepare_eqtl_for_plotting,
    validate_eqtl_df,
)

# Validation utilities
from .utils import ValidationError, to_pandas

__all__ = [
    # Core
    "__version__",
    "LocusZoomPlotter",
    # Backends
    "BackendType",
    "get_backend",
    # Reference data
    "download_dog_recombination_maps",
    # Colors
    "get_ld_color",
    "get_ld_bin",
    "get_ld_color_palette",
    "LEAD_SNP_COLOR",
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
    # Logging
    "enable_logging",
    "disable_logging",
    # Validation & Utils
    "ValidationError",
    "to_pandas",
]

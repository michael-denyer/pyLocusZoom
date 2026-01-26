"""pyLocusZoom - Regional association plots for GWAS results.

This package provides LocusZoom-style regional association plots with:
- LD coloring based on R² with lead variant
- Gene and exon tracks
- Recombination rate overlays (dog built-in, or user-provided)
- Automatic SNP labeling

Example:
    >>> from pylocuszoom import LocusZoomPlotter
    >>> plotter = LocusZoomPlotter(species="dog")
    >>> fig = plotter.plot(gwas_df, chrom=1, start=1000000, end=2000000)
    >>> fig.savefig("regional_plot.png", dpi=150)

Species Support:
    - Dog: Full features including built-in recombination maps
    - Cat: LD coloring and gene tracks (user provides recombination data)
    - Custom: User provides all reference data
"""

__version__ = "0.1.0"

# Main plotter class
# Utility functions that may be useful standalone
from .colors import LEAD_SNP_COLOR, get_ld_bin, get_ld_color, get_ld_color_palette
from .gene_track import get_nearest_gene, plot_gene_track
from .labels import add_snp_labels
from .ld import calculate_ld

# Logging configuration
from .logging import disable_logging, enable_logging
from .plotter import LocusZoomPlotter

# Reference data management
from .recombination import (
    add_recombination_overlay,
    download_dog_recombination_maps,
    get_recombination_rate_for_region,
    load_recombination_map,
)

# Validation utilities
from .utils import ValidationError

__all__ = [
    # Core
    "__version__",
    "LocusZoomPlotter",
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
    # Logging
    "enable_logging",
    "disable_logging",
    # Validation
    "ValidationError",
]

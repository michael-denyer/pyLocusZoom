"""Miami plot generator for comparing two GWAS datasets.

Provides visualization of GWAS comparisons with mirrored y-axes:
- Top panel shows -log10(p) ascending (standard Manhattan)
- Bottom panel shows -log10(p) descending (inverted y-axis)
- Both panels share x-axis with consistent chromosome alignment
"""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._family_renderers import MiamiRenderer
from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
)
from .backends import BackendType, get_backend
from .manhattan import prepare_manhattan_data


class MiamiPlotter:
    """Miami plot generator for comparing two GWAS datasets.

    Creates mirrored Manhattan plots with top panel showing -log10(p)
    ascending and bottom panel showing -log10(p) descending, enabling
    visual comparison of two GWAS results.

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        species: Species name ('canine', 'feline', 'human', or None).
            Used to determine chromosome order.
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').
        genomewide_threshold: P-value threshold for significance line.

    Example:
        >>> plotter = MiamiPlotter(species="human")
        >>> fig = plotter.plot_miami(discovery_df, replication_df)
        >>> fig.savefig("miami_plot.png", dpi=150)
    """

    def __init__(
        self,
        species: str = "canine",
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
    ):
        """Initialize the Miami plotter."""
        self.species = species
        self._backend = get_backend(backend)
        self._renderer = MiamiRenderer(self._backend)
        self.backend_name = backend
        self.genomewide_threshold = genomewide_threshold

    def plot_miami(
        self,
        top_df: pd.DataFrame,
        bottom_df: pd.DataFrame,
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        rs_col: Optional[str] = None,
        custom_chrom_order: Optional[List[str]] = None,
        top_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        bottom_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        top_label: Optional[str] = None,
        bottom_label: Optional[str] = None,
        top_snp_annotations: Optional[List[str]] = None,
        bottom_snp_annotations: Optional[List[str]] = None,
        highlight_regions: Optional[List[Tuple[str, int, int]]] = None,
        highlight_color: str = "yellow",
        highlight_alpha: float = 0.3,
        figsize: Tuple[float, float] = (12, 8),
        title: Optional[str] = None,
    ) -> Any:
        """Create a Miami plot comparing two GWAS datasets.

        The top panel displays -log10(p) values ascending (standard Manhattan),
        while the bottom panel displays -log10(p) descending (inverted), creating
        a mirrored comparison.

        Args:
            top_df: GWAS results DataFrame for top panel.
            bottom_df: GWAS results DataFrame for bottom panel.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP RS ID (for hover tooltips and annotations).
            custom_chrom_order: Custom chromosome order (overrides species).
            top_threshold: Significance threshold for top panel. None to skip.
            bottom_threshold: Significance threshold for bottom panel. None to skip.
            top_label: Label for top panel (e.g., "Discovery").
            bottom_label: Label for bottom panel (e.g., "Replication").
            top_snp_annotations: List of SNP IDs to annotate on top panel.
                Requires rs_col to be set. Basic text labels (no collision avoidance).
            bottom_snp_annotations: List of SNP IDs to annotate on bottom panel.
                Requires rs_col to be set. Basic text labels (no collision avoidance).
            highlight_regions: List of (chrom, start, end) tuples to highlight.
                Regions are drawn as vertical spans across both panels.
            highlight_color: Color for highlighted regions.
            highlight_alpha: Transparency for highlighted regions (0-1).
            figsize: Figure size as (width, height).
            title: Overall plot title.

        Returns:
            Figure object (type depends on backend).

        Raises:
            ValueError: If required columns are missing from either DataFrame.

        Example:
            >>> fig = plotter.plot_miami(
            ...     discovery_df,
            ...     replication_df,
            ...     top_label="Discovery",
            ...     bottom_label="Replication",
            ... )
        """
        # Compute union of chromosomes to ensure consistent alignment
        # This is critical to avoid Pitfall #3 from research
        all_chroms = self._get_chromosome_union(top_df, bottom_df, chrom_col)

        # Use custom order if provided, otherwise use chromosome union
        if custom_chrom_order is None:
            custom_chrom_order = all_chroms

        # Prepare both datasets with consistent chromosome ordering
        top_prepared = prepare_manhattan_data(
            df=top_df,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )
        bottom_prepared = prepare_manhattan_data(
            df=bottom_df,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )

        return self._renderer.render(
            top_prepared,
            bottom_prepared,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
            top_threshold=top_threshold,
            bottom_threshold=bottom_threshold,
            top_label=top_label,
            bottom_label=bottom_label,
            top_snp_annotations=top_snp_annotations,
            bottom_snp_annotations=bottom_snp_annotations,
            highlight_regions=highlight_regions,
            highlight_color=highlight_color,
            highlight_alpha=highlight_alpha,
            figsize=figsize,
            title=title,
        )

    def _get_chromosome_union(
        self,
        top_df: pd.DataFrame,
        bottom_df: pd.DataFrame,
        chrom_col: str,
    ) -> List[str]:
        """Get union of chromosomes from both DataFrames.

        Ensures consistent chromosome ordering across both panels,
        which is critical for x-axis alignment in Miami plots.

        Args:
            top_df: Top panel DataFrame.
            bottom_df: Bottom panel DataFrame.
            chrom_col: Chromosome column name.

        Returns:
            Sorted list of all unique chromosomes.
        """
        top_chroms = set(top_df[chrom_col].astype(str).unique())
        bottom_chroms = set(bottom_df[chrom_col].astype(str).unique())
        all_chroms = top_chroms | bottom_chroms

        # Sort chromosomes: numeric first (by value), then alphabetic
        def sort_key(chrom: str) -> tuple:
            try:
                return (0, int(chrom), "")
            except ValueError:
                return (1, 0, chrom)

        return sorted(all_chroms, key=sort_key)

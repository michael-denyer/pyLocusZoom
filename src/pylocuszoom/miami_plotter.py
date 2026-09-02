"""Miami plot generator for comparing two GWAS datasets.

Provides visualization of GWAS comparisons with mirrored y-axes:
- Top panel shows -log10(p) ascending (standard Manhattan)
- Bottom panel shows -log10(p) descending (inverted y-axis)
- Both panels share x-axis with consistent chromosome alignment
"""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._figure import render_figure
from ._miami_panels import MiamiRequest, miami_plan
from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from .backends import BackendType, get_backend
from .backends.hover import HoverConfig
from .manhattan import prepare_manhattan_frames
from .species import Species, resolve_species


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
        species: Species name, alias or record ('canine', 'dog', 'feline',
            'human', or None). An unknown name raises ValidationError.
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
        species: str | Species | None = "canine",
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
    ):
        """Initialize the Miami plotter."""
        self.species = resolve_species(species)
        self._backend = get_backend(backend)
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
        top_threshold: ThresholdArg = UNSET,
        bottom_threshold: ThresholdArg = UNSET,
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
                Both panels are laid out against it.
            top_threshold: Significance threshold for the top panel. Defaults to
                the plotter's ``genomewide_threshold``; pass None to draw no line.
            bottom_threshold: Significance threshold for the bottom panel. Same
                defaulting as ``top_threshold``.
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
            ValueError: If required columns are missing from either DataFrame,
                or if the plotter's species is unknown.

        Example:
            >>> fig = plotter.plot_miami(
            ...     discovery_df,
            ...     replication_df,
            ...     top_label="Discovery",
            ...     bottom_label="Replication",
            ... )
        """
        top_threshold = resolve_threshold(top_threshold, self.genomewide_threshold)
        bottom_threshold = resolve_threshold(
            bottom_threshold, self.genomewide_threshold
        )

        top_prepared, bottom_prepared = prepare_manhattan_frames(
            [top_df, bottom_df],
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )

        request = MiamiRequest(
            top=top_prepared,
            bottom=bottom_prepared,
            hover=(
                HoverConfig(snp_col=rs_col, pos_col=pos_col, p_col=p_col)
                if rs_col is not None
                else None
            ),
            rs_col=rs_col,
            top_threshold=top_threshold,
            bottom_threshold=bottom_threshold,
            top_label=top_label,
            bottom_label=bottom_label,
            top_annotations=tuple(top_snp_annotations or ()),
            bottom_annotations=tuple(bottom_snp_annotations or ()),
            highlights=tuple(highlight_regions or ()),
            highlight_color=highlight_color,
            highlight_alpha=highlight_alpha,
            figsize=figsize,
            title=title,
        )
        return render_figure(self._backend, miami_plan(request))

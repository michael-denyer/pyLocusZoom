"""Miami plot generator for comparing two GWAS datasets.

Provides visualization of GWAS comparisons with mirrored y-axes:
- Top panel shows -log10(p) ascending (standard Manhattan)
- Bottom panel shows -log10(p) descending (inverted y-axis)
- Both panels share x-axis with consistent chromosome alignment
"""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
    add_significance_line,
)
from .backends import BackendType, get_backend
from .backends.hover import HoverConfig, HoverDataBuilder
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
            rs_col: Column name for SNP RS ID (for hover tooltips).
            custom_chrom_order: Custom chromosome order (overrides species).
            top_threshold: Significance threshold for top panel.
            bottom_threshold: Significance threshold for bottom panel.
            top_label: Label for top panel (e.g., "Discovery").
            bottom_label: Label for bottom panel (e.g., "Replication").
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

        # Create figure with 2 panels
        fig, axes = self._backend.create_figure(
            n_panels=2,
            height_ratios=[1.0, 1.0],
            figsize=figsize,
            sharex=True,
        )
        top_ax, bottom_ax = axes

        # Get consistent chrom order from first prepared dataframe
        chrom_order = top_prepared.attrs["chrom_order"]
        chrom_centers = top_prepared.attrs["chrom_centers"]

        # Plot top panel (normal y-axis)
        self._render_manhattan_points(
            ax=top_ax,
            prepared_df=top_prepared,
            chrom_order=chrom_order,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
        )
        add_significance_line(self._backend, top_ax, top_threshold)

        # Plot bottom panel (inverted y-axis)
        self._render_manhattan_points(
            ax=bottom_ax,
            prepared_df=bottom_prepared,
            chrom_order=chrom_order,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
        )
        add_significance_line(self._backend, bottom_ax, bottom_threshold)

        # Set x limits consistently across both panels
        x_min = min(
            top_prepared["_cumulative_pos"].min(),
            bottom_prepared["_cumulative_pos"].min(),
        )
        x_max = max(
            top_prepared["_cumulative_pos"].max(),
            bottom_prepared["_cumulative_pos"].max(),
        )
        x_padding = (x_max - x_min) * 0.01
        self._backend.set_xlim(top_ax, x_min - x_padding, x_max + x_padding)
        self._backend.set_xlim(bottom_ax, x_min - x_padding, x_max + x_padding)

        # Set y limits - critical for Miami plot
        top_y_max = top_prepared["_neg_log_p"].max() * 1.1
        bottom_y_max = bottom_prepared["_neg_log_p"].max() * 1.1

        # Top panel: normal y-axis (0 at bottom, max at top)
        self._backend.set_ylim(top_ax, 0, top_y_max)

        # Bottom panel: inverted y-axis (max at bottom, 0 at top)
        # For matplotlib, passing (max, 0) inverts the axis
        self._backend.set_ylim(bottom_ax, bottom_y_max, 0)

        # Set x-axis ticks for chromosome labels
        positions = [
            chrom_centers[chrom] for chrom in chrom_order if chrom in chrom_centers
        ]
        labels = [str(chrom) for chrom in chrom_order if chrom in chrom_centers]

        # Both panels get x-axis ticks (needed for interactive backends)
        self._backend.set_xticks(top_ax, positions, labels, fontsize=8)
        self._backend.set_xticks(bottom_ax, positions, labels, fontsize=8)

        # Y-axis labels
        self._backend.set_ylabel(top_ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_ylabel(bottom_ax, r"$-\log_{10}(p)$", fontsize=12)

        # X-axis label only on bottom panel
        self._backend.set_xlabel(bottom_ax, "Chromosome", fontsize=12)

        # Panel labels
        if top_label:
            self._backend.add_panel_label(top_ax, top_label)
        if bottom_label:
            self._backend.add_panel_label(bottom_ax, bottom_label)

        # Hide spines for clean appearance
        self._backend.hide_spines(top_ax, ["top", "right"])
        self._backend.hide_spines(bottom_ax, ["top", "right"])

        # Overall title
        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.92, hspace=0.05)
        else:
            self._backend.finalize_layout(fig, hspace=0.05)

        return fig

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

    def _render_manhattan_points(
        self,
        ax: Any,
        prepared_df: pd.DataFrame,
        chrom_order: List[str],
        chrom_col: str,
        pos_col: str,
        p_col: str,
        rs_col: Optional[str] = None,
        point_size: int = MANHATTAN_POINT_SIZE,
    ) -> None:
        """Render Manhattan plot scatter points grouped by chromosome.

        Args:
            ax: Axes object from backend.
            prepared_df: DataFrame with _chrom_str, _cumulative_pos, _neg_log_p, _color.
            chrom_order: List of chromosome names in display order.
            chrom_col: Original chromosome column name.
            pos_col: Original position column name.
            p_col: Original p-value column name.
            rs_col: RS ID column name for hover data.
            point_size: Size of scatter points.
        """
        for chrom in chrom_order:
            chrom_data = prepared_df[prepared_df["_chrom_str"] == chrom]
            if len(chrom_data) > 0:
                # Build hover data for interactive backends
                hover_df = None
                if self._backend.supports_hover and rs_col is not None:
                    hover_config = HoverConfig(
                        snp_col=rs_col,
                        pos_col=pos_col,
                        p_col=p_col,
                    )
                    builder = HoverDataBuilder(hover_config)
                    hover_df = builder.build_dataframe(chrom_data)

                self._backend.scatter(
                    ax,
                    chrom_data["_cumulative_pos"],
                    chrom_data["_neg_log_p"],
                    colors=chrom_data["_color"].iloc[0],
                    sizes=point_size,
                    marker="o",
                    edgecolor=POINT_EDGE_COLOR,
                    linewidth=MANHATTAN_EDGE_WIDTH,
                    zorder=2,
                    hover_data=hover_df,
                )

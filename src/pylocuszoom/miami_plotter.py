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


# [1c:MiamiPlotter] Two-trait mirrored Manhattan plots — see docs/CODEMAP.md
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

        # Set y limits - critical for Miami plot.
        # Guard against degenerate ylim(0, 0) when all p=1 or DataFrames empty.
        raw_top = top_prepared["_neg_log_p"].max()
        raw_bot = bottom_prepared["_neg_log_p"].max()
        top_y_max = max(raw_top * 1.1, 1.0) if pd.notna(raw_top) else 1.0
        bottom_y_max = max(raw_bot * 1.1, 1.0) if pd.notna(raw_bot) else 1.0

        # Top panel: normal y-axis (0 at bottom, max at top)
        self._backend.set_ylim(top_ax, 0, top_y_max)

        # Bottom panel: inverted y-axis (max at bottom, 0 at top)
        # For matplotlib, passing (max, 0) inverts the axis
        self._backend.set_ylim(bottom_ax, bottom_y_max, 0)

        # Set x-axis ticks for chromosome labels (needed for both panels in interactive backends)
        valid_chroms = [c for c in chrom_order if c in chrom_centers]
        positions = [chrom_centers[c] for c in valid_chroms]
        labels = [str(c) for c in valid_chroms]
        for ax in [top_ax, bottom_ax]:
            self._backend.set_xticks(ax, positions, labels, fontsize=8)

        # Y-axis labels
        self._backend.set_ylabel(top_ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_ylabel(bottom_ax, r"$-\log_{10}(p)$", fontsize=12)

        # X-axis label only on bottom panel
        self._backend.set_xlabel(bottom_ax, "Chromosome", fontsize=12)

        # Panel labels - top at top, bottom at bottom for Miami plot layout
        if top_label:
            self._backend.add_panel_label(top_ax, top_label, y_frac=0.95)
        if bottom_label:
            # For Miami plots, bottom panel label should be at the bottom of the panel
            self._backend.add_panel_label(bottom_ax, bottom_label, y_frac=0.05)

        # SNP annotations
        if top_snp_annotations and rs_col:
            self._add_snp_annotations(
                ax=top_ax,
                prepared_df=top_prepared,
                rs_col=rs_col,
                snp_ids=top_snp_annotations,
            )
        if bottom_snp_annotations and rs_col:
            self._add_snp_annotations(
                ax=bottom_ax,
                prepared_df=bottom_prepared,
                rs_col=rs_col,
                snp_ids=bottom_snp_annotations,
            )

        # Region highlighting
        if highlight_regions:
            # Calculate chromosome offsets for position conversion
            chrom_offsets = self._get_chrom_offsets(top_prepared, pos_col)
            for chrom, start, end in highlight_regions:
                self._draw_region_highlight(
                    fig=fig,
                    top_ax=top_ax,
                    bottom_ax=bottom_ax,
                    chrom=str(chrom),
                    start=start,
                    end=end,
                    chrom_offsets=chrom_offsets,
                    color=highlight_color,
                    alpha=highlight_alpha,
                )

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

    def _add_snp_annotations(
        self,
        ax: Any,
        prepared_df: pd.DataFrame,
        rs_col: str,
        snp_ids: List[str],
    ) -> None:
        """Add text annotations for specified SNPs.

        Basic annotation without collision avoidance. For matplotlib,
        use add_snp_labels() in LocusZoomPlotter for adjustText support.

        Args:
            ax: Axes object from backend.
            prepared_df: DataFrame with _cumulative_pos, _neg_log_p, and RS column.
            rs_col: Column name containing SNP IDs.
            snp_ids: List of SNP IDs to annotate.
        """
        # Filter to only requested SNPs
        snps_to_annotate = prepared_df[prepared_df[rs_col].isin(snp_ids)]

        for _, row in snps_to_annotate.iterrows():
            x = row["_cumulative_pos"]
            y = row["_neg_log_p"]
            label = row[rs_col]

            # Add text slightly above the point
            self._backend.add_text(
                ax,
                x=x,
                y=y,
                text=str(label),
                fontsize=8,
                ha="center",
                va="bottom",
            )

    def _get_chrom_offsets(
        self, prepared_df: pd.DataFrame, pos_col: str
    ) -> dict[str, float]:
        """Calculate cumulative position offset for each chromosome.

        The offset is the difference between cumulative position and original
        position for the first SNP on each chromosome.

        Args:
            prepared_df: DataFrame with _chrom_str, _cumulative_pos columns.
            pos_col: Original position column name.

        Returns:
            Dict mapping chromosome string to its cumulative offset.
        """
        offsets = {}
        for chrom in prepared_df.attrs.get("chrom_order", []):
            chrom_data = prepared_df[prepared_df["_chrom_str"] == str(chrom)]
            if not chrom_data.empty:
                first_row = chrom_data.iloc[0]
                offsets[str(chrom)] = first_row["_cumulative_pos"] - first_row[pos_col]
        return offsets

    def _draw_region_highlight(
        self,
        fig: Any,
        top_ax: Any,
        bottom_ax: Any,
        chrom: str,
        start: int,
        end: int,
        chrom_offsets: dict[str, float],
        color: str,
        alpha: float,
    ) -> None:
        """Draw highlighted region across both panels.

        Uses backend-specific implementations since region highlighting
        is not part of the PlotBackend protocol.

        Args:
            fig: Figure object from backend.
            top_ax: Top panel axes.
            bottom_ax: Bottom panel axes.
            chrom: Chromosome name (as string).
            start: Region start position (bp).
            end: Region end position (bp).
            chrom_offsets: Dict mapping chromosome to cumulative offset.
            color: Highlight color.
            alpha: Highlight transparency.
        """
        if chrom not in chrom_offsets:
            return  # Chromosome not in data

        offset = chrom_offsets[chrom]
        x_start = offset + start
        x_end = offset + end

        if self.backend_name == "matplotlib":
            top_ax.axvspan(x_start, x_end, color=color, alpha=alpha, zorder=0)
            bottom_ax.axvspan(x_start, x_end, color=color, alpha=alpha, zorder=0)

        elif self.backend_name == "plotly":
            # For plotly, fig is the Figure and axes are (fig, row) tuples
            fig.add_vrect(
                x0=x_start,
                x1=x_end,
                fillcolor=color,
                opacity=alpha,
                layer="below",
                line_width=0,
                row=1,
                col=1,
            )
            fig.add_vrect(
                x0=x_start,
                x1=x_end,
                fillcolor=color,
                opacity=alpha,
                layer="below",
                line_width=0,
                row=2,
                col=1,
            )

        elif self.backend_name == "bokeh":
            from bokeh.models import BoxAnnotation

            # For bokeh, axes are figure objects
            for ax in [top_ax, bottom_ax]:
                box = BoxAnnotation(
                    left=x_start,
                    right=x_end,
                    fill_color=color,
                    fill_alpha=alpha,
                )
                ax.add_layout(box)

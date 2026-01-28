"""Base protocol for plotting backends.

Defines the interface that matplotlib, plotly, and bokeh backends must implement.
"""

from typing import Any, Callable, List, Optional, Protocol, Tuple, Union

import pandas as pd


class PlotBackend(Protocol):
    """Protocol defining the backend interface for LocusZoom plots.

    All backends (matplotlib, plotly, bokeh) must implement these methods
    to enable consistent plotting across different rendering engines.

    Capability Properties:
        supports_snp_labels: Whether backend supports text labels via adjustText.
        supports_hover: Whether backend supports hover tooltips.
        supports_secondary_axis: Whether backend supports twin y-axis for overlays.
    """

    # =========================================================================
    # Capability Properties
    # =========================================================================

    @property
    def supports_snp_labels(self) -> bool:
        """Whether backend supports text labels via adjustText.

        Matplotlib supports SNP labels using adjustText for automatic repositioning.
        Interactive backends (Plotly, Bokeh) use hover tooltips instead.
        """
        ...

    @property
    def supports_hover(self) -> bool:
        """Whether backend supports hover tooltips.

        Interactive backends (Plotly, Bokeh) support hover tooltips.
        Matplotlib does not support hover - use SNP labels instead.
        """
        ...

    @property
    def supports_secondary_axis(self) -> bool:
        """Whether backend supports twin y-axis for recombination overlay.

        All current backends support secondary axes, but this allows for
        future backends that may not.
        """
        ...

    # =========================================================================
    # Figure Creation
    # =========================================================================

    def create_figure(
        self,
        n_panels: int,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[Any, List[Any]]:
        """Create a figure with multiple panels (subplots).

        Args:
            n_panels: Number of vertical panels.
            height_ratios: Relative heights for each panel.
            figsize: Figure size as (width, height).
            sharex: Whether panels share the x-axis.

        Returns:
            Tuple of (figure, list of axes/panels).
        """
        ...

    def finalize_layout(
        self,
        fig: Any,
        left: float = 0.08,
        right: float = 0.95,
        top: float = 0.95,
        bottom: float = 0.1,
        hspace: float = 0.08,
    ) -> None:
        """Finalize figure layout with margins and spacing.

        Args:
            fig: Figure object.
            left: Left margin fraction.
            right: Right margin fraction.
            top: Top margin fraction.
            bottom: Bottom margin fraction.
            hspace: Vertical space between subplots.
        """
        ...

    # =========================================================================
    # Basic Plotting
    # =========================================================================

    def scatter(
        self,
        ax: Any,
        x: pd.Series,
        y: pd.Series,
        colors: Union[str, List[str], pd.Series],
        sizes: Union[float, List[float], pd.Series] = 60,
        marker: str = "o",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
        hover_data: Optional[pd.DataFrame] = None,
        label: Optional[str] = None,
    ) -> Any:
        """Create a scatter plot on the given axes.

        Args:
            ax: Axes or panel to plot on.
            x: X-axis values (positions).
            y: Y-axis values (-log10 p-values).
            colors: Point colors (single color or per-point).
            sizes: Point sizes.
            marker: Marker style.
            edgecolor: Marker edge color.
            linewidth: Marker edge width.
            zorder: Drawing order.
            hover_data: DataFrame with columns for hover tooltips.
            label: Legend label.

        Returns:
            The scatter plot object.
        """
        ...

    def line(
        self,
        ax: Any,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        zorder: int = 1,
        label: Optional[str] = None,
    ) -> Any:
        """Create a line plot on the given axes.

        Args:
            ax: Axes or panel to plot on.
            x: X-axis values.
            y: Y-axis values.
            color: Line color.
            linewidth: Line width.
            alpha: Transparency.
            linestyle: Line style ('-', '--', ':', '-.').
            zorder: Drawing order.
            label: Legend label.

        Returns:
            The line plot object.
        """
        ...

    def fill_between(
        self,
        ax: Any,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> Any:
        """Fill area between two y-values.

        Args:
            ax: Axes or panel to plot on.
            x: X-axis values.
            y1: Lower y boundary.
            y2: Upper y boundary.
            color: Fill color.
            alpha: Transparency.
            zorder: Drawing order.

        Returns:
            The fill object.
        """
        ...

    def axhline(
        self,
        ax: Any,
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the axes.

        Args:
            ax: Axes or panel.
            y: Y-value for the line.
            color: Line color.
            linestyle: Line style.
            linewidth: Line width.
            alpha: Line transparency (0-1).
            zorder: Drawing order.

        Returns:
            The line object.
        """
        ...

    def axvline(
        self,
        ax: Any,
        x: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a vertical line across the axes.

        Args:
            ax: Axes or panel.
            x: X-value for the line.
            color: Line color.
            linestyle: Line style.
            linewidth: Line width.
            alpha: Line transparency (0-1).
            zorder: Drawing order.

        Returns:
            The line object.
        """
        ...

    # =========================================================================
    # Text and Annotations
    # =========================================================================

    def add_text(
        self,
        ax: Any,
        x: float,
        y: float,
        text: str,
        fontsize: int = 10,
        ha: str = "center",
        va: str = "bottom",
        rotation: float = 0,
        color: str = "black",
    ) -> Any:
        """Add text annotation to axes.

        Args:
            ax: Axes or panel.
            x: X position.
            y: Y position.
            text: Text content.
            fontsize: Font size.
            ha: Horizontal alignment.
            va: Vertical alignment.
            rotation: Text rotation in degrees.
            color: Text color.

        Returns:
            The text object.
        """
        ...

    def add_panel_label(
        self,
        ax: Any,
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel.

        Used for panel letters (A, B, C) in multi-panel figures.

        Args:
            ax: Axes or panel.
            label: Label text (e.g., "A", "B").
            x_frac: Horizontal position as fraction of axes (0-1).
            y_frac: Vertical position as fraction of axes (0-1).
        """
        ...

    def add_snp_labels(
        self,
        ax: Any,
        df: pd.DataFrame,
        pos_col: str,
        neglog10p_col: str,
        rs_col: str,
        label_top_n: int,
        genes_df: Optional[pd.DataFrame],
        chrom: int,
    ) -> None:
        """Add SNP labels to plot.

        No-op if supports_snp_labels=False. Matplotlib uses adjustText
        for automatic label repositioning to avoid overlaps.

        Args:
            ax: Axes or panel.
            df: DataFrame with SNP data.
            pos_col: Column name for position.
            neglog10p_col: Column name for -log10(p-value).
            rs_col: Column name for SNP ID.
            label_top_n: Number of top SNPs to label.
            genes_df: Gene annotations (unused, for signature compatibility).
            chrom: Chromosome number (unused, for signature compatibility).
        """
        ...

    # =========================================================================
    # Shapes and Patches
    # =========================================================================

    def add_rectangle(
        self,
        ax: Any,
        xy: Tuple[float, float],
        width: float,
        height: float,
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a rectangle patch to axes.

        Args:
            ax: Axes or panel.
            xy: Bottom-left corner coordinates.
            width: Rectangle width.
            height: Rectangle height.
            facecolor: Fill color.
            edgecolor: Edge color.
            linewidth: Edge width.
            zorder: Drawing order.

        Returns:
            The rectangle object.
        """
        ...

    def add_polygon(
        self,
        ax: Any,
        points: List[List[float]],
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add polygon patch to axes.

        Used for gene track directional arrows.

        Args:
            ax: Axes or panel.
            points: List of [x, y] coordinate pairs forming the polygon.
            facecolor: Fill color.
            edgecolor: Edge color.
            linewidth: Edge width.
            zorder: Drawing order.

        Returns:
            The polygon object.
        """
        ...

    # =========================================================================
    # Axis Configuration
    # =========================================================================

    def set_xlim(self, ax: Any, left: float, right: float) -> None:
        """Set x-axis limits.

        Args:
            ax: Axes or panel.
            left: Minimum x value.
            right: Maximum x value.
        """
        ...

    def set_ylim(self, ax: Any, bottom: float, top: float) -> None:
        """Set y-axis limits.

        Args:
            ax: Axes or panel.
            bottom: Minimum y value.
            top: Maximum y value.
        """
        ...

    def set_xlabel(self, ax: Any, label: str, fontsize: int = 12) -> None:
        """Set x-axis label.

        Args:
            ax: Axes or panel.
            label: Label text.
            fontsize: Font size.
        """
        ...

    def set_ylabel(self, ax: Any, label: str, fontsize: int = 12) -> None:
        """Set y-axis label.

        Args:
            ax: Axes or panel.
            label: Label text.
            fontsize: Font size.
        """
        ...

    def set_yticks(
        self,
        ax: Any,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
    ) -> None:
        """Set y-axis tick positions and labels.

        Args:
            ax: Axes or panel.
            positions: Tick positions.
            labels: Tick labels.
            fontsize: Font size.
        """
        ...

    def set_title(self, ax: Any, title: str, fontsize: int = 14) -> None:
        """Set panel title.

        Args:
            ax: Axes or panel.
            title: Title text.
            fontsize: Font size.
        """
        ...

    def hide_spines(self, ax: Any, spines: List[str]) -> None:
        """Hide specified axis spines.

        Args:
            ax: Axes or panel.
            spines: List of spine names ('top', 'right', 'bottom', 'left').
        """
        ...

    def hide_yaxis(self, ax: Any) -> None:
        """Hide y-axis for gene track panels.

        Hides y-axis ticks, labels, and line. Gene tracks don't need
        a y-axis since the vertical position is just for layout.

        Args:
            ax: Axes or panel.
        """
        ...

    def format_xaxis_mb(self, ax: Any) -> None:
        """Format x-axis to show megabase values.

        Args:
            ax: Axes or panel.
        """
        ...

    # =========================================================================
    # Secondary Y-Axis (for recombination overlay)
    # =========================================================================

    def create_twin_axis(self, ax: Any) -> Any:
        """Create a secondary y-axis sharing the same x-axis.

        Args:
            ax: Primary axes.

        Returns:
            Secondary axes for overlay (e.g., recombination rate).
        """
        ...

    def line_secondary(
        self,
        ax: Any,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        label: Optional[str] = None,
        yaxis_name: Any = None,
    ) -> Any:
        """Create line on secondary y-axis.

        Args:
            ax: Axes or panel (may be tuple for Plotly).
            x: X-axis values.
            y: Y-axis values.
            color: Line color.
            linewidth: Line width.
            alpha: Transparency.
            linestyle: Line style.
            label: Legend label.
            yaxis_name: Backend-specific secondary axis identifier.

        Returns:
            The line object.
        """
        ...

    def fill_between_secondary(
        self,
        ax: Any,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        yaxis_name: Any = None,
    ) -> Any:
        """Fill area on secondary y-axis.

        Args:
            ax: Axes or panel.
            x: X-axis values.
            y1: Lower y boundary.
            y2: Upper y boundary.
            color: Fill color.
            alpha: Transparency.
            yaxis_name: Backend-specific secondary axis identifier.

        Returns:
            The fill object.
        """
        ...

    def set_secondary_ylim(
        self,
        ax: Any,
        bottom: float,
        top: float,
        yaxis_name: Any = None,
    ) -> None:
        """Set secondary y-axis limits.

        Args:
            ax: Axes or panel.
            bottom: Minimum y value.
            top: Maximum y value.
            yaxis_name: Backend-specific secondary axis identifier.
        """
        ...

    def set_secondary_ylabel(
        self,
        ax: Any,
        label: str,
        color: str = "black",
        fontsize: int = 10,
        yaxis_name: Any = None,
    ) -> None:
        """Set secondary y-axis label.

        Args:
            ax: Axes or panel.
            label: Label text.
            color: Label color.
            fontsize: Font size.
            yaxis_name: Backend-specific secondary axis identifier.
        """
        ...

    # =========================================================================
    # Legends
    # =========================================================================

    def add_legend(
        self,
        ax: Any,
        handles: List[Any],
        labels: List[str],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> Any:
        """Add a legend to the axes.

        Args:
            ax: Axes or panel.
            handles: Legend handle objects.
            labels: Legend labels.
            loc: Legend location.
            title: Legend title.

        Returns:
            The legend object.
        """
        ...

    def add_ld_legend(
        self,
        ax: Any,
        ld_bins: List[Tuple[float, str, str]],
        lead_snp_color: str,
    ) -> None:
        """Add LD color legend.

        Shows the linkage disequilibrium (r^2) color scale and lead SNP marker.

        Args:
            ax: Axes or panel.
            ld_bins: List of (threshold, label, color) tuples defining LD bins.
            lead_snp_color: Color for lead SNP marker in legend.
        """
        ...

    def add_eqtl_legend(
        self,
        ax: Any,
        eqtl_positive_bins: List[Tuple[float, float, str, str]],
        eqtl_negative_bins: List[Tuple[float, float, str, str]],
    ) -> None:
        """Add eQTL effect size legend to the axes.

        Args:
            ax: Axes or panel.
            eqtl_positive_bins: List of (min, max, label, color) for positive effects.
            eqtl_negative_bins: List of (min, max, label, color) for negative effects.
        """
        ...

    def add_finemapping_legend(
        self,
        ax: Any,
        credible_sets: List[int],
        get_color_func: Callable[[int], str],
    ) -> None:
        """Add fine-mapping credible set legend to the axes.

        Args:
            ax: Axes or panel.
            credible_sets: List of credible set IDs to show.
            get_color_func: Function that takes CS ID and returns color.
        """
        ...

    def add_simple_legend(
        self,
        ax: Any,
        label: str,
        loc: str = "upper right",
    ) -> None:
        """Add a simple legend entry for labeled data already in the plot.

        Args:
            ax: Axes or panel.
            label: Legend label for labeled scatter data.
            loc: Legend location.
        """
        ...

    # =========================================================================
    # Specialized Charts
    # =========================================================================

    def hbar(
        self,
        ax: Any,
        y: pd.Series,
        width: pd.Series,
        height: float = 0.8,
        left: Union[float, pd.Series] = 0,
        color: Union[str, List[str]] = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Create horizontal bar chart.

        Args:
            ax: Axes or panel.
            y: Y positions for bars.
            width: Bar widths (x-extent).
            height: Bar height.
            left: Left edge positions.
            color: Bar colors.
            edgecolor: Edge color.
            linewidth: Edge width.
            zorder: Drawing order.

        Returns:
            The bar collection object.
        """
        ...

    def errorbar_h(
        self,
        ax: Any,
        x: pd.Series,
        y: pd.Series,
        xerr_lower: pd.Series,
        xerr_upper: pd.Series,
        color: str = "black",
        linewidth: float = 1.5,
        capsize: float = 3,
        zorder: int = 3,
    ) -> Any:
        """Add horizontal error bars (for forest plots).

        Args:
            ax: Axes or panel.
            x: X positions (effect sizes).
            y: Y positions.
            xerr_lower: Lower error (distance from x).
            xerr_upper: Upper error (distance from x).
            color: Line color.
            linewidth: Line width.
            capsize: Cap size in points.
            zorder: Drawing order.

        Returns:
            The errorbar object.
        """
        ...

    # =========================================================================
    # File Operations
    # =========================================================================

    def save(
        self,
        fig: Any,
        path: str,
        dpi: int = 150,
        bbox_inches: str = "tight",
    ) -> None:
        """Save figure to file.

        Args:
            fig: Figure object.
            path: Output file path (.png, .pdf, .html).
            dpi: Resolution for raster formats.
            bbox_inches: Bounding box adjustment.
        """
        ...

    def show(self, fig: Any) -> None:
        """Display the figure.

        Args:
            fig: Figure object.
        """
        ...

    def close(self, fig: Any) -> None:
        """Close the figure and free resources.

        Args:
            fig: Figure object.
        """
        ...

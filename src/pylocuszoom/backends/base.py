"""Base protocol for plotting backends.

Defines the interface that matplotlib, plotly, and bokeh backends must implement.
"""

from typing import Any, List, Optional, Protocol, Tuple, Union

import pandas as pd


class PlotBackend(Protocol):
    """Protocol defining the backend interface for LocusZoom plots.

    All backends (matplotlib, plotly, bokeh) must implement these methods
    to enable consistent plotting across different rendering engines.
    """

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

    def set_title(self, ax: Any, title: str, fontsize: int = 14) -> None:
        """Set panel title.

        Args:
            ax: Axes or panel.
            title: Title text.
            fontsize: Font size.
        """
        ...

    def create_twin_axis(self, ax: Any) -> Any:
        """Create a secondary y-axis sharing the same x-axis.

        Args:
            ax: Primary axes.

        Returns:
            Secondary axes for overlay (e.g., recombination rate).
        """
        ...

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

    def hide_spines(self, ax: Any, spines: List[str]) -> None:
        """Hide specified axis spines.

        Args:
            ax: Axes or panel.
            spines: List of spine names ('top', 'right', 'bottom', 'left').
        """
        ...

    def format_xaxis_mb(self, ax: Any) -> None:
        """Format x-axis to show megabase values.

        Args:
            ax: Axes or panel.
        """
        ...

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

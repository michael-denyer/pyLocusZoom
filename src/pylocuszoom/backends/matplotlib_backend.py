"""Matplotlib backend for pyLocusZoom.

Default backend providing static publication-quality plots.
"""

from typing import Any, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator

from . import register_backend
from .composition import LegendEntry, heatmap_highlight_cells


@register_backend("matplotlib")
class MatplotlibBackend:
    """Matplotlib backend for static plot generation.

    This is the default backend, producing publication-quality static plots
    suitable for papers and presentations.

    Capability Properties:
        supports_snp_labels: True - uses adjustText for automatic label positioning.
        supports_hover: False - static plots don't support hover tooltips.
        supports_secondary_axis: True - supports twin y-axis via twinx().
    """

    # =========================================================================
    # Capability Properties
    # =========================================================================

    @property
    def supports_hover(self) -> bool:
        """Matplotlib does not support hover tooltips."""
        return False

    def create_figure(
        self,
        n_panels: int,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[Figure, List[Axes]]:
        """Create a figure with multiple panels."""
        if n_panels == 1:
            fig, ax = plt.subplots(figsize=figsize)
            self._hide_top_right(ax)
            return fig, [ax]

        fig, axes = plt.subplots(
            n_panels,
            1,
            figsize=figsize,
            height_ratios=height_ratios,
            sharex=sharex,
        )

        for ax in axes:
            self._hide_top_right(ax)
        return fig, list(axes)

    def create_figure_grid(
        self,
        n_rows: int,
        n_cols: int,
        width_ratios: Optional[List[float]] = None,
        height_ratios: Optional[List[float]] = None,
        figsize: Tuple[float, float] = (12.0, 8.0),
    ) -> Tuple[Figure, List[Axes]]:
        """Create a figure with a grid of subplots."""
        gridspec_kw = {}
        if width_ratios is not None:
            gridspec_kw["width_ratios"] = width_ratios
        if height_ratios is not None:
            gridspec_kw["height_ratios"] = height_ratios

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=figsize,
            gridspec_kw=gridspec_kw if gridspec_kw else None,
        )

        # Flatten axes to list
        import numpy as np

        flat = list(axes.flatten()) if isinstance(axes, np.ndarray) else [axes]
        for ax in flat:
            self._hide_top_right(ax)
        return fig, flat

    @staticmethod
    def _hide_top_right(ax: Axes) -> None:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def scatter(
        self,
        ax: Axes,
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

        Note: hover_data is ignored for matplotlib (static plots).
        """
        return ax.scatter(
            x,
            y,
            c=colors,
            s=sizes,
            marker=marker,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
            label=label,
        )

    def line(
        self,
        ax: Axes,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        zorder: int = 1,
        label: Optional[str] = None,
    ) -> Any:
        """Create a line plot on the given axes."""
        (line,) = ax.plot(
            x,
            y,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            linestyle=linestyle,
            zorder=zorder,
            label=label,
        )
        return line

    def fill_between(
        self,
        ax: Axes,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> Any:
        """Fill area between two y-values."""
        return ax.fill_between(x, y1, y2, color=color, alpha=alpha, zorder=zorder)

    def axhline(
        self,
        ax: Axes,
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the axes."""
        return ax.axhline(
            y=y,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
        )

    def add_text(
        self,
        ax: Axes,
        x: float,
        y: float,
        text: str,
        fontsize: int = 10,
        ha: str = "center",
        va: str = "bottom",
        rotation: float = 0,
        color: str = "black",
    ) -> Any:
        """Add text annotation to axes."""
        return ax.text(
            x, y, text, fontsize=fontsize, ha=ha, va=va, rotation=rotation, color=color
        )

    def add_panel_label(
        self,
        ax: Axes,
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        ax.annotate(
            label,
            xy=(x_frac, y_frac),
            xycoords="axes fraction",
            fontsize=10,
            fontweight="bold",
            ha="left",
            va="top",
        )

    def add_snp_labels(
        self,
        ax: Axes,
        df: pd.DataFrame,
        pos_col: str,
        neglog10p_col: str,
        rs_col: str,
        label_top_n: int,
        adjust: bool = True,
        lead_pos: Optional[int] = None,
        region_span: Optional[int] = None,
    ) -> List[Any]:
        """Add SNP labels using adjustText."""
        from ..labels import add_snp_labels as _add_snp_labels

        return _add_snp_labels(
            ax,
            df,
            pos_col=pos_col,
            neglog10p_col=neglog10p_col,
            rs_col=rs_col,
            label_top_n=label_top_n,
            adjust=adjust,
            lead_pos=lead_pos,
            region_span=region_span,
        )

    def add_rectangle(
        self,
        ax: Axes,
        xy: Tuple[float, float],
        width: float,
        height: float,
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a rectangle patch to axes."""
        rect = Rectangle(
            xy,
            width,
            height,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
        ax.add_patch(rect)
        return rect

    def add_polygon(
        self,
        ax: Axes,
        points: List[List[float]],
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a polygon patch to axes."""
        polygon = Polygon(
            points,
            closed=True,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
        ax.add_patch(polygon)
        return polygon

    def set_xlim(self, ax: Axes, left: float, right: float) -> None:
        """Set x-axis limits."""
        ax.set_xlim(left, right)

    def set_ylim(self, ax: Axes, bottom: float, top: float) -> None:
        """Set y-axis limits."""
        ax.set_ylim(bottom, top)

    def set_xlabel(self, ax: Axes, label: str, fontsize: int = 12) -> None:
        """Set x-axis label."""
        ax.set_xlabel(label, fontsize=fontsize)

    def set_ylabel(self, ax: Axes, label: str, fontsize: int = 12) -> None:
        """Set y-axis label."""
        ax.set_ylabel(label, fontsize=fontsize)

    def set_yticks(
        self,
        ax: Axes,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
    ) -> None:
        """Set y-axis tick positions and labels."""
        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=fontsize)

    def set_xticks(
        self,
        ax: Axes,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
        rotation: int = 0,
        ha: str = "center",
    ) -> None:
        """Set x-axis tick positions and labels."""
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=fontsize, rotation=rotation, ha=ha)

    def set_title(self, ax: Axes, title: str, fontsize: int = 14) -> None:
        """Set panel title."""
        ax.set_title(
            title,
            fontsize=fontsize,
            fontweight="bold",
            fontfamily="sans-serif",
        )

    def set_suptitle(self, fig: Figure, title: str, fontsize: int = 14) -> None:
        """Set overall figure title (super title)."""
        fig.suptitle(title, fontsize=fontsize, fontweight="bold")

    def create_twin_axis(self, ax: Axes) -> Axes:
        """Create a secondary y-axis sharing the same x-axis."""
        secondary = ax.twinx()
        ax.spines["right"].set_visible(True)
        secondary.spines["top"].set_visible(False)
        return secondary

    def line_secondary(
        self,
        secondary: Axes,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        label: Optional[str] = None,
    ) -> Any:
        """Create a line on the secondary (twin) axes."""
        return self.line(
            secondary,
            x,
            y,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            linestyle=linestyle,
            label=label,
        )

    def fill_between_secondary(
        self,
        secondary: Axes,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
    ) -> Any:
        """Fill area on the secondary (twin) axes."""
        return self.fill_between(secondary, x, y1, y2, color=color, alpha=alpha)

    def set_secondary_ylim(
        self,
        secondary: Axes,
        bottom: float,
        top: float,
    ) -> None:
        """Set secondary y-axis limits."""
        self.set_ylim(secondary, bottom, top)

    def set_secondary_ylabel(
        self,
        secondary: Axes,
        label: str,
        color: str = "black",
        fontsize: int = 10,
    ) -> None:
        """Set secondary y-axis label."""
        secondary.set_ylabel(label, fontsize=fontsize, color=color)
        secondary.tick_params(axis="y", labelcolor=color, labelsize=fontsize - 1)

    def add_legend(
        self,
        ax: Axes,
        entries: List[LegendEntry],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> Any:
        """Render backend-neutral legend entries as matplotlib handles."""
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        handles: List[Any] = []
        for entry in entries:
            edge = entry.edgecolor or "black"
            if entry.marker == "patch":
                handles.append(
                    Patch(facecolor=entry.color, edgecolor=edge, label=entry.label)
                )
            else:
                handles.append(
                    Line2D(
                        [0],
                        [0],
                        marker=entry.marker,
                        color="w",
                        markerfacecolor=entry.color,
                        markeredgecolor=edge,
                        markersize=7,
                        label=entry.label,
                    )
                )
        return ax.legend(
            handles=handles,
            loc=loc,
            title=title,
            fontsize=9,
            frameon=True,
            framealpha=0.9,
            title_fontsize=10,
            handlelength=1.5,
            handleheight=1.0,
            labelspacing=0.4,
        )

    def hide_yaxis(self, ax: Axes) -> None:
        """Hide y-axis ticks, labels, and line."""
        ax.yaxis.set_visible(False)
        ax.spines["left"].set_visible(False)

    def format_xaxis_mb(self, ax: Axes) -> None:
        """Format x-axis to show megabase values."""
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:.2f}"))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))

    def axvline(
        self,
        ax: Axes,
        x: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a vertical line across the axes."""
        return ax.axvline(
            x=x,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
        )

    def hbar(
        self,
        ax: Axes,
        y: pd.Series,
        width: pd.Series,
        height: float = 0.8,
        left: Union[float, pd.Series] = 0,
        color: Union[str, List[str]] = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Create horizontal bar chart."""
        return ax.barh(
            y=y,
            width=width,
            height=height,
            left=left,
            color=color,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )

    def errorbar_h(
        self,
        ax: Axes,
        x: pd.Series,
        y: pd.Series,
        xerr_lower: pd.Series,
        xerr_upper: pd.Series,
        color: str = "black",
        linewidth: float = 1.5,
        capsize: float = 3,
        zorder: int = 3,
    ) -> Any:
        """Add horizontal error bars."""
        xerr = [xerr_lower.values, xerr_upper.values]
        return ax.errorbar(
            x=x,
            y=y,
            xerr=xerr,
            fmt="none",
            ecolor=color,
            elinewidth=linewidth,
            capsize=capsize,
            zorder=zorder,
        )

    def finalize_layout(
        self,
        fig: Figure,
        left: float = 0.08,
        right: float = 0.95,
        top: float = 0.95,
        bottom: float = 0.1,
        hspace: float = 0.08,
    ) -> None:
        """Adjust subplot layout parameters."""
        fig.subplots_adjust(
            left=left, right=right, top=top, bottom=bottom, hspace=hspace
        )

    def add_region_highlight(
        self,
        fig: Figure,
        axes: List[Axes],
        x_start: float,
        x_end: float,
        color: str = "yellow",
        alpha: float = 0.3,
    ) -> None:
        """Highlight an x-range across multiple matplotlib axes."""
        for ax in axes:
            ax.axvspan(x_start, x_end, color=color, alpha=alpha, zorder=0)

    def highlight_heatmap_snp(
        self,
        ax: Axes,
        fig: Figure,
        snp_idx: int,
        n_snps: int,
        color: str = "#FF0000",
        linewidth: float = 2,
    ) -> None:
        """Highlight a SNP's row and column in the heatmap."""
        for x, y in heatmap_highlight_cells(snp_idx, n_snps):
            ax.add_patch(
                Rectangle(
                    (x - 0.5, y - 0.5),
                    1.0,
                    1.0,
                    fill=False,
                    edgecolor=color,
                    linewidth=linewidth,
                    zorder=10,
                )
            )

    def add_heatmap(
        self,
        ax: Axes,
        data: Any,
        x_coords: List[float],
        y_coords: List[float],
        cmap_colors: Optional[List[str]] = None,
        vmin: float = 0.0,
        vmax: float = 1.0,
    ) -> Any:
        """Render a heatmap of an already-shaped matrix."""
        from matplotlib.colors import LinearSegmentedColormap

        # Default white-to-red gradient
        if cmap_colors is None:
            cmap_colors = ["#FFFFFF", "#FF0000"]

        cmap = LinearSegmentedColormap.from_list("ld_heatmap", cmap_colors, N=256)

        # Compute extent from coordinates if non-trivial (genomic coordinates)
        # This allows heatmap to align with regional plot x-axis
        extent = None
        if x_coords and len(x_coords) > 1:
            x_min, x_max = min(x_coords), max(x_coords)
            # Only use extent if coordinates are not simple indices
            if x_max - x_min > len(x_coords):
                # Genomic coordinates - use extent for alignment
                # Add half-cell padding for proper cell centering
                y_min, y_max = min(y_coords), max(y_coords)
                extent = [x_min, x_max, y_min - 0.5, y_max + 0.5]

        # Use imshow for heatmap rendering
        im = ax.imshow(
            data,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
            origin="lower",
            extent=extent,
        )
        return im

    def add_colorbar(
        self,
        ax: Axes,
        mappable: Any,
        label: str = "R²",
        orientation: str = "vertical",
    ) -> Any:
        """Add colorbar legend for heatmap."""
        cbar = plt.colorbar(mappable, ax=ax, orientation=orientation)
        cbar.set_label(label)
        return cbar

"""Matplotlib backend for pyLocusZoom.

Default backend providing static publication-quality plots.
"""

from typing import Any, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator

from ..colors import FOOTER_COLOR
from . import register_backend
from .composition import LegendEntry, cell_edges

# Side and bottom margins, as fractions of the figure. No caller has ever
# varied them, so they are this backend's own layout policy rather than part
# of the finalize_layout contract.
_LEFT_MARGIN = 0.08
_RIGHT_MARGIN = 0.95
_BOTTOM_MARGIN = 0.1
_INSET_COLORBAR_LABEL = "_pylocuszoom_colorbar"


@register_backend("matplotlib")
class MatplotlibBackend:
    """Matplotlib backend for static plot generation.

    This is the default backend, producing publication-quality static plots
    suitable for papers and presentations. It is the one built-in backend
    that implements ``SupportsSNPLabels``, through adjustText.
    """

    def create_figure(
        self,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[Figure, List[Axes]]:
        """Create a figure with one panel per height ratio."""
        if len(height_ratios) == 1:
            fig, ax = plt.subplots(figsize=figsize)
            self._hide_top_right(ax)
            return fig, [ax]

        fig, axes = plt.subplots(
            len(height_ratios),
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
        alpha: Optional[float] = None,
    ) -> None:
        """Create a scatter plot on the given axes.

        Note: hover_data is ignored for matplotlib (static plots).
        """
        ax.scatter(
            x,
            y,
            c=colors,
            s=sizes,
            marker=marker,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
            alpha=alpha,
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
    ) -> None:
        """Create a line plot on the given axes."""
        ax.plot(
            x,
            y,
            color=color,
            linewidth=linewidth,
            alpha=alpha,
            linestyle=linestyle,
            zorder=zorder,
        )

    def fill_between(
        self,
        ax: Axes,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> None:
        """Fill area between two y-values."""
        ax.fill_between(x, y1, y2, color=color, alpha=alpha, zorder=zorder)

    def axhline(
        self,
        ax: Axes,
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> None:
        """Add a horizontal line across the axes."""
        ax.axhline(
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
        color: str = "black",
    ) -> None:
        """Add text annotation to axes."""
        ax.text(x, y, text, fontsize=fontsize, ha=ha, va=va, color=color)

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
    ) -> None:
        """Add SNP labels using adjustText."""
        from ..labels import add_snp_labels as _add_snp_labels

        _add_snp_labels(
            ax,
            df,
            pos_col=pos_col,
            neglog10p_col=neglog10p_col,
            rs_col=rs_col,
            label_top_n=label_top_n,
        )

    def add_rectangle(
        self,
        ax: Axes,
        xy: Tuple[float, float],
        width: float,
        height: float,
        facecolor: Optional[str] = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> None:
        """Add a rectangle patch to axes."""
        rect = Rectangle(
            xy,
            width,
            height,
            fill=facecolor is not None,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
        ax.add_patch(rect)

    def add_polygon(
        self,
        ax: Axes,
        points: List[List[float]],
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> None:
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

    def set_tick_fontsize(self, ax: Axes, fontsize: int) -> None:
        """Set the tick label size on both axes."""
        ax.tick_params(axis="both", labelsize=fontsize)

    def set_title(
        self,
        ax: Axes,
        title: str,
        fontsize: int = 14,
        fontweight: Literal["bold", "normal"] = "bold",
    ) -> None:
        """Set panel title."""
        ax.set_title(
            title,
            fontsize=fontsize,
            fontweight=fontweight,
            fontfamily="sans-serif",
        )

    def set_suptitle(
        self,
        fig: Figure,
        title: str,
        fontsize: int = 14,
        fontweight: Literal["bold", "normal"] = "bold",
    ) -> None:
        """Set overall figure title (super title)."""
        fig.suptitle(title, fontsize=fontsize, fontweight=fontweight)

    def set_footer(self, fig: Figure, text: str, fontsize: int = 10) -> None:
        """Raise the panels by one text line and write the footer beneath."""
        line = 2 * fontsize / 72 / fig.get_figheight()
        fig.subplots_adjust(bottom=fig.subplotpars.bottom + line)
        fig.text(
            0.5,
            0.005,
            text,
            ha="center",
            va="bottom",
            fontsize=fontsize,
            style="italic",
            color=FOOTER_COLOR,
            usetex=False,
            parse_math=False,
        )

    def create_twin_axis(self, ax: Axes) -> Axes:
        """Create a secondary y-axis sharing the same x-axis."""
        secondary = ax.twinx()
        ax.spines["right"].set_visible(True)
        secondary.spines["top"].set_visible(False)
        return secondary

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
        title: Optional[str] = None,
    ) -> None:
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
        ax.legend(
            handles=handles,
            loc="upper right",
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
    ) -> None:
        """Add a vertical line across the axes."""
        ax.axvline(
            x=x,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
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
    ) -> None:
        """Add horizontal error bars."""
        xerr = [xerr_lower.values, xerr_upper.values]
        ax.errorbar(
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
        top: float = 0.95,
        hspace: float = 0.08,
    ) -> None:
        """Adjust subplot layout parameters."""
        # A colorbar outside a shared-x panel needs figure-wide space so the
        # gene, association and heatmap axes retain identical genomic scales.
        has_inset_colorbar = any(
            ax.get_label() == _INSET_COLORBAR_LABEL for ax in fig.axes
        )
        fig.subplots_adjust(
            left=_LEFT_MARGIN,
            right=0.88 if has_inset_colorbar else _RIGHT_MARGIN,
            top=top,
            bottom=_BOTTOM_MARGIN,
            hspace=hspace,
        )

    def add_region_highlight(
        self,
        axes: List[Axes],
        x_start: float,
        x_end: float,
        color: str = "yellow",
        alpha: float = 0.3,
    ) -> None:
        """Highlight an x-range across multiple matplotlib axes."""
        for ax in axes:
            ax.axvspan(x_start, x_end, color=color, alpha=alpha, zorder=0)

    def add_heatmap(
        self,
        ax: Axes,
        data: Any,
        x_coords: List[float],
        y_coords: List[float],
        cmap_colors: List[str],
        vmin: float = 0.0,
        vmax: float = 1.0,
        colorbar_label: Optional[str] = None,
    ) -> None:
        """Render a heatmap of an already-shaped matrix."""
        from matplotlib.colors import LinearSegmentedColormap

        cmap = LinearSegmentedColormap.from_list("ld_heatmap", cmap_colors, N=256)

        x_edges = cell_edges(x_coords)
        y_edges = cell_edges(y_coords)
        mesh = ax.pcolormesh(
            [x_edges[0][0], *(right for _, right in x_edges)],
            [y_edges[0][0], *(right for _, right in y_edges)],
            data,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            shading="flat",
            rasterized=True,
        )
        if colorbar_label is None:
            return
        if len(ax.get_shared_x_axes().get_siblings(ax)) > 1:
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes

            # Native colorbar allocation shrinks only its parent subplot. An
            # inset follows that subplot after finalize_layout, without changing
            # its width relative to the other genomic panels.
            colorbar_ax = inset_axes(
                ax,
                width="2%",
                height="100%",
                loc="lower left",
                bbox_to_anchor=(1.02, 0, 1, 1),
                bbox_transform=ax.transAxes,
                borderpad=0,
            )
            cbar = ax.figure.colorbar(mesh, cax=colorbar_ax)
            colorbar_ax.set_label(_INSET_COLORBAR_LABEL)
        else:
            cbar = ax.figure.colorbar(mesh, ax=ax)
        cbar.set_label(colorbar_label)

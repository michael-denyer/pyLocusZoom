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
    def supports_snp_labels(self) -> bool:
        """Matplotlib supports SNP labels via adjustText."""
        return True

    @property
    def supports_hover(self) -> bool:
        """Matplotlib does not support hover tooltips."""
        return False

    @property
    def supports_secondary_axis(self) -> bool:
        """Matplotlib supports twin y-axis."""
        return True

    def create_figure(
        self,
        n_panels: int,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[Figure, List[Axes]]:
        """Create a figure with multiple panels.

        Args:
            n_panels: Number of vertical panels.
            height_ratios: Relative heights for each panel.
            figsize: Figure size as (width, height).
            sharex: Whether panels share the x-axis.

        Returns:
            Tuple of (figure, list of axes).
        """
        # Prevent auto-display in interactive environments
        plt.ioff()

        if n_panels == 1:
            fig, ax = plt.subplots(figsize=figsize)
            return fig, [ax]

        fig, axes = plt.subplots(
            n_panels,
            1,
            figsize=figsize,
            height_ratios=height_ratios,
            sharex=sharex,
        )

        return fig, list(axes)

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
        """Add label text at fractional position in panel.

        Args:
            ax: Matplotlib axes.
            label: Label text (e.g., "A", "B").
            x_frac: Horizontal position as fraction of axes (0-1).
            y_frac: Vertical position as fraction of axes (0-1).
        """
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
        genes_df: Optional[pd.DataFrame],
        chrom: int,
    ) -> None:
        """Add SNP labels using adjustText.

        Args:
            ax: Matplotlib axes.
            df: DataFrame with SNP data.
            pos_col: Column name for position.
            neglog10p_col: Column name for -log10(p-value).
            rs_col: Column name for SNP ID.
            label_top_n: Number of top SNPs to label.
            genes_df: Gene annotations (unused, for signature compatibility).
            chrom: Chromosome number (unused, for signature compatibility).
        """
        from ..labels import add_snp_labels as _add_snp_labels

        _add_snp_labels(
            ax,
            df,
            pos_col=pos_col,
            neglog10p_col=neglog10p_col,
            rs_col=rs_col,
            label_top_n=label_top_n,
            genes_df=genes_df,
            chrom=chrom,
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

    def create_twin_axis(self, ax: Axes) -> Axes:
        """Create a secondary y-axis sharing the same x-axis."""
        return ax.twinx()

    def line_secondary(
        self,
        ax: Axes,
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

        For matplotlib, the ax should already be a twin axis from create_twin_axis().
        The yaxis_name parameter is ignored (provided for interface compatibility).

        Args:
            ax: Secondary axes from create_twin_axis().
            x: X-axis values.
            y: Y-axis values.
            color: Line color.
            linewidth: Line width.
            alpha: Transparency.
            linestyle: Line style.
            label: Legend label.
            yaxis_name: Ignored for matplotlib.

        Returns:
            The line object.
        """
        return self.line(
            ax,
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
        ax: Axes,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        yaxis_name: Any = None,
    ) -> Any:
        """Fill area on secondary y-axis.

        For matplotlib, the ax should already be a twin axis from create_twin_axis().
        The yaxis_name parameter is ignored (provided for interface compatibility).

        Args:
            ax: Secondary axes from create_twin_axis().
            x: X-axis values.
            y1: Lower y boundary.
            y2: Upper y boundary.
            color: Fill color.
            alpha: Transparency.
            yaxis_name: Ignored for matplotlib.

        Returns:
            The fill object.
        """
        return self.fill_between(ax, x, y1, y2, color=color, alpha=alpha)

    def set_secondary_ylim(
        self,
        ax: Axes,
        bottom: float,
        top: float,
        yaxis_name: Any = None,
    ) -> None:
        """Set secondary y-axis limits.

        For matplotlib, the ax should already be a twin axis from create_twin_axis().
        The yaxis_name parameter is ignored (provided for interface compatibility).

        Args:
            ax: Secondary axes from create_twin_axis().
            bottom: Minimum y value.
            top: Maximum y value.
            yaxis_name: Ignored for matplotlib.
        """
        self.set_ylim(ax, bottom, top)

    def set_secondary_ylabel(
        self,
        ax: Axes,
        label: str,
        color: str = "black",
        fontsize: int = 10,
        yaxis_name: Any = None,
    ) -> None:
        """Set secondary y-axis label.

        For matplotlib, the ax should already be a twin axis from create_twin_axis().
        The yaxis_name parameter is ignored (provided for interface compatibility).

        Args:
            ax: Secondary axes from create_twin_axis().
            label: Label text.
            color: Label color.
            fontsize: Font size.
            yaxis_name: Ignored for matplotlib.
        """
        ax.set_ylabel(label, fontsize=fontsize, color=color)
        ax.tick_params(axis="y", labelcolor=color, labelsize=fontsize - 1)

    def add_legend(
        self,
        ax: Axes,
        handles: List[Any],
        labels: List[str],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> Any:
        """Add a legend to the axes."""
        return ax.legend(
            handles=handles,
            labels=labels,
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

    def hide_spines(self, ax: Axes, spines: List[str]) -> None:
        """Hide specified axis spines."""
        for spine in spines:
            ax.spines[spine].set_visible(False)

    def hide_yaxis(self, ax: Axes) -> None:
        """Hide y-axis ticks, labels, and line."""
        ax.yaxis.set_visible(False)

    def format_xaxis_mb(self, ax: Axes) -> None:
        """Format x-axis to show megabase values."""
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:.2f}"))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))

    def save(
        self,
        fig: Figure,
        path: str,
        dpi: int = 150,
        bbox_inches: str = "tight",
    ) -> None:
        """Save figure to file."""
        fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches)

    def show(self, fig: Figure) -> None:
        """Display the figure."""
        plt.ion()
        plt.show()

    def close(self, fig: Figure) -> None:
        """Close the figure and free resources."""
        plt.close(fig)

    def add_eqtl_legend(
        self,
        ax: Axes,
        eqtl_positive_bins: List[Tuple[float, float, str, str]],
        eqtl_negative_bins: List[Tuple[float, float, str, str]],
    ) -> None:
        """Add eQTL effect size legend using matplotlib Line2D markers."""
        from matplotlib.lines import Line2D

        legend_elements = []

        # Positive effects (upward triangles)
        for _, _, label, color in eqtl_positive_bins:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="^",
                    color="w",
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markersize=7,
                    label=label,
                )
            )

        # Negative effects (downward triangles)
        for _, _, label, color in eqtl_negative_bins:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="v",
                    color="w",
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markersize=7,
                    label=label,
                )
            )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=8,
            frameon=True,
            framealpha=0.9,
            title="eQTL effect",
            title_fontsize=9,
            handlelength=1.2,
            handleheight=1.0,
            labelspacing=0.3,
        )

    def add_finemapping_legend(
        self,
        ax: Axes,
        credible_sets: List[int],
        get_color_func: Any,
    ) -> None:
        """Add fine-mapping credible set legend using matplotlib Line2D markers."""
        from matplotlib.lines import Line2D

        if not credible_sets:
            return

        legend_elements = []
        for cs_id in credible_sets:
            color = get_color_func(cs_id)
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markersize=7,
                    label=f"CS{cs_id}",
                )
            )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=8,
            frameon=True,
            framealpha=0.9,
            title="Credible sets",
            title_fontsize=9,
            handlelength=1.2,
            handleheight=1.0,
            labelspacing=0.3,
        )

    def add_simple_legend(
        self,
        ax: Axes,
        label: str,
        loc: str = "upper right",
    ) -> None:
        """Add simple legend for labeled scatter data."""
        ax.legend(loc=loc, fontsize=9)

    def add_ld_legend(
        self,
        ax: Axes,
        ld_bins: List[Tuple[float, str, str]],
        lead_snp_color: str,
    ) -> None:
        """Add LD color legend using matplotlib patches.

        Args:
            ax: Matplotlib axes.
            ld_bins: List of (threshold, label, color) tuples defining LD bins.
            lead_snp_color: Color for lead SNP marker in legend.
        """
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        from ..colors import get_ld_color_palette

        palette = get_ld_color_palette()
        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="D",
                color="w",
                markerfacecolor=lead_snp_color,
                markeredgecolor="black",
                markersize=6,
                label="Lead SNP",
            ),
        ]
        for _threshold, label, _color in ld_bins:
            legend_elements.append(
                Patch(facecolor=palette[label], edgecolor="black", label=label)
            )
        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=9,
            frameon=True,
            framealpha=0.9,
            title=r"$r^2$",
            title_fontsize=10,
            handlelength=1.5,
            handleheight=1.0,
            labelspacing=0.4,
        )

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
        """Adjust subplot layout parameters.

        Args:
            fig: Figure object.
            left: Left margin.
            right: Right margin.
            top: Top margin.
            bottom: Bottom margin.
            hspace: Height space between subplots.
        """
        fig.subplots_adjust(
            left=left, right=right, top=top, bottom=bottom, hspace=hspace
        )
        plt.ion()

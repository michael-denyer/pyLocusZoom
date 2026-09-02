"""Plotly backend for pyLocusZoom.

Interactive backend with hover tooltips and zoom/pan capabilities.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from . import convert_latex_to_unicode, register_backend
from ._coerce import (
    broadcast,
    marker_colors,
    marker_diameter,
    normalize_ratios,
    pixels,
)
from .composition import LegendEntry, heatmap_highlight_cells, mb_tick_positions
from .hover import plotly_hovertemplate
from .plotly_layout import (
    _Panel,
    configure_legend,
    secondary_axis_key,
    x_range,
)

# Style mappings (matplotlib -> Plotly)
_MARKER_SYMBOLS = {
    "o": "circle",
    "D": "diamond",
    "s": "square",
    "^": "triangle-up",
    "v": "triangle-down",
}
_DASH_MAP = {
    "-": "solid",
    "--": "dash",
    ":": "dot",
    "-.": "dashdot",
}


@register_backend("plotly")
class PlotlyBackend:
    """Plotly backend for interactive plot generation.

    Produces interactive HTML plots with hover tooltips showing:
    - SNP RS ID
    - P-value
    - R² with lead SNP
    - Nearest gene
    """

    @property
    def supports_hover(self) -> bool:
        """Plotly supports hover tooltips."""
        return True

    def create_figure(
        self,
        n_panels: int,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[go.Figure, List[Any]]:
        """Create a figure with multiple panels."""
        width_px, height_px = pixels(figsize)
        row_heights = normalize_ratios(height_ratios)

        fig = make_subplots(
            rows=n_panels,
            cols=1,
            shared_xaxes=sharex,
            vertical_spacing=0.02,
            row_heights=row_heights,
        )

        fig.update_layout(
            width=width_px,
            height=height_px,
            showlegend=True,
            template="plotly_white",
        )

        # Style all panels for clean LocusZoom appearance
        axis_style = dict(
            showgrid=False,
            showline=True,
            linecolor="black",
            ticks="outside",
            minor_ticks="",
            zeroline=False,
        )
        for row in range(1, n_panels + 1):
            panel = _Panel(fig, row)
            fig.update_layout(
                **{panel.axis("xaxis"): axis_style, panel.axis("yaxis"): axis_style}
            )

        # Return (fig, row) tuples for each panel
        # This matches the expected ax parameter format for all methods
        panel_refs = [(fig, row) for row in range(1, n_panels + 1)]
        return fig, panel_refs

    def create_figure_grid(
        self,
        n_rows: int,
        n_cols: int,
        width_ratios: Optional[List[float]] = None,
        height_ratios: Optional[List[float]] = None,
        figsize: Tuple[float, float] = (12.0, 8.0),
    ) -> Tuple[go.Figure, List[Any]]:
        """Create a figure with a grid of subplots."""
        width_px, height_px = pixels(figsize)
        column_widths = normalize_ratios(width_ratios)
        row_heights_norm = normalize_ratios(height_ratios)

        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            column_widths=column_widths,
            row_heights=row_heights_norm,
            horizontal_spacing=0.08,
            vertical_spacing=0.08,
        )

        fig.update_layout(
            width=width_px,
            height=height_px,
            showlegend=True,
            template="plotly_white",
        )

        # Style all panels
        axis_style = dict(
            showgrid=False,
            showline=True,
            linecolor="black",
            ticks="outside",
            zeroline=False,
        )
        for row in range(1, n_rows + 1):
            for col in range(1, n_cols + 1):
                panel = _Panel(fig, row, col, n_cols)
                fig.update_layout(
                    **{panel.axis("xaxis"): axis_style, panel.axis("yaxis"): axis_style}
                )

        # Return flattened list of (fig, row, col, n_cols) tuples
        panel_refs = []
        for row in range(1, n_rows + 1):
            for col in range(1, n_cols + 1):
                panel_refs.append((fig, row, col, n_cols))
        return fig, panel_refs

    def scatter(
        self,
        ax: Tuple[go.Figure, int],
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
        """Create a scatter plot on the given panel.

        For plotly, ax is a tuple of (figure, row_number) or (figure, row, col).
        """
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        # Convert matplotlib marker to plotly symbol
        symbol = _MARKER_SYMBOLS.get(marker, "circle")

        size = marker_diameter(sizes)

        # Build hover template
        if hover_data is not None:
            customdata = hover_data.values
            hovertemplate = plotly_hovertemplate(hover_data)
        else:
            customdata = None
            hovertemplate = "x: %{x}<br>y: %{y:.2f}<extra></extra>"

        marker_color = marker_colors(colors)

        trace = go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(
                color=marker_color,
                size=size,
                symbol=symbol,
                line=dict(color=edgecolor, width=linewidth),
            ),
            customdata=customdata,
            hovertemplate=hovertemplate,
            name=label or "",
            showlegend=label is not None,
        )

        fig.add_trace(trace, row=row, col=col)
        return trace

    def line(
        self,
        ax: Tuple[go.Figure, int],
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        zorder: int = 1,
        label: Optional[str] = None,
    ) -> Any:
        """Create a line plot on the given panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col
        dash = _DASH_MAP.get(linestyle, "solid")

        trace = go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color=color, width=linewidth, dash=dash),
            opacity=alpha,
            name=label or "",
            showlegend=label is not None,
        )

        fig.add_trace(trace, row=row, col=col)
        return trace

    def fill_between(
        self,
        ax: Tuple[go.Figure, int],
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> Any:
        """Fill area between two y-values."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        y1 = pd.Series(broadcast(y1, len(x)))

        trace = go.Scatter(
            x=pd.concat([x, x[::-1]]),
            y=pd.concat([y2, y1[::-1]]),
            fill="toself",
            fillcolor=color,
            opacity=alpha,
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
        )

        fig.add_trace(trace, row=row, col=col)
        return trace

    def axhline(
        self,
        ax: Tuple[go.Figure, int],
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col
        dash = _DASH_MAP.get(linestyle, "dash")

        fig.add_hline(
            y=y,
            line_dash=dash,
            line_color=color,
            line_width=linewidth,
            opacity=alpha,
            row=row,
            col=col,
        )

    def add_text(
        self,
        ax: Tuple[go.Figure, int],
        x: float,
        y: float,
        text: str,
        fontsize: int = 10,
        ha: str = "center",
        va: str = "bottom",
        rotation: float = 0,
        color: str = "black",
    ) -> Any:
        """Add text annotation to panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        # Map alignment
        xanchor_map = {"center": "center", "left": "left", "right": "right"}
        yanchor_map = {"bottom": "bottom", "top": "top", "center": "middle"}

        fig.add_annotation(
            x=x,
            y=y,
            text=text,
            font=dict(size=fontsize, color=color),
            xanchor=xanchor_map.get(ha, "center"),
            yanchor=yanchor_map.get(va, "bottom"),
            textangle=-rotation,
            showarrow=False,
            row=row,
            col=col,
        )

    def add_rectangle(
        self,
        ax: Tuple[go.Figure, int],
        xy: Tuple[float, float],
        width: float,
        height: float,
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a rectangle to the panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        x0, y0 = xy
        x1, y1 = x0 + width, y0 + height

        fig.add_shape(
            type="rect",
            x0=x0,
            y0=y0,
            x1=x1,
            y1=y1,
            fillcolor=facecolor,
            line=dict(color=edgecolor, width=linewidth),
            row=row,
            col=col,
        )

    def add_polygon(
        self,
        ax: Tuple[go.Figure, int],
        points: List[List[float]],
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a polygon (e.g., triangle for strand arrows) to the panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        # Build SVG path from points
        path = f"M {points[0][0]} {points[0][1]}"
        for px, py in points[1:]:
            path += f" L {px} {py}"
        path += " Z"

        fig.add_shape(
            type="path",
            path=path,
            fillcolor=facecolor,
            line=dict(color=edgecolor, width=linewidth),
            row=row,
            col=col,
        )

    def set_xlim(self, ax: Tuple[go.Figure, int], left: float, right: float) -> None:
        """Set x-axis limits."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        fig.update_layout(**{panel.axis("xaxis"): dict(range=[left, right])})

    def set_ylim(self, ax: Tuple[go.Figure, int], bottom: float, top: float) -> None:
        """Set y-axis limits."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        fig.update_layout(**{panel.axis("yaxis"): dict(range=[bottom, top])})

    def set_xlabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set x-axis label."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        label = convert_latex_to_unicode(label)
        fig.update_layout(
            **{
                panel.axis("xaxis"): dict(
                    title=dict(text=label, font=dict(size=fontsize))
                )
            }
        )

    def set_ylabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set y-axis label."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        label = convert_latex_to_unicode(label)
        fig.update_layout(
            **{
                panel.axis("yaxis"): dict(
                    title=dict(text=label, font=dict(size=fontsize))
                )
            }
        )

    def set_yticks(
        self,
        ax: Tuple[go.Figure, int],
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
    ) -> None:
        """Set y-axis tick positions and labels."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        fig.update_layout(
            **{
                panel.axis("yaxis"): dict(
                    tickmode="array",
                    tickvals=positions,
                    ticktext=labels,
                    tickfont=dict(size=fontsize),
                )
            }
        )

    def set_xticks(
        self,
        ax: Tuple[go.Figure, int],
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
        rotation: int = 0,
        ha: str = "center",
    ) -> None:
        """Set x-axis tick positions and labels."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        fig.update_layout(
            **{
                panel.axis("xaxis"): dict(
                    tickmode="array",
                    tickvals=positions,
                    ticktext=labels,
                    tickfont=dict(size=fontsize),
                    tickangle=-rotation if rotation else 0,
                )
            }
        )

    def set_title(
        self, ax: Tuple[go.Figure, int], title: str, fontsize: int = 14
    ) -> None:
        """Set subplot title using annotation.

        For grid layouts, this adds an annotation above the subplot.
        For single-column layouts, sets the global figure title for the first panel.
        """
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel

        if n_cols == 1 and row == 1:
            # Single-column layout: use global figure title
            fig.update_layout(title=dict(text=title, font=dict(size=fontsize)))
        else:
            # Grid layout: add annotation above the subplot
            # Use subplot's axis domain for positioning
            xref = f"{panel.ref('x')} domain"
            yref = f"{panel.ref('y')} domain"

            fig.add_annotation(
                text=f"<b>{title}</b>",
                xref=xref,
                yref=yref,
                x=0.5,
                y=1.05,
                showarrow=False,
                font=dict(size=fontsize),
                xanchor="center",
                yanchor="bottom",
            )

    def set_suptitle(self, fig: go.Figure, title: str, fontsize: int = 14) -> None:
        """Set overall figure title (super title)."""
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=fontsize),
                x=0.5,
                xanchor="center",
            )
        )

    def create_twin_axis(self, ax: Tuple[go.Figure, int]) -> Any:
        """Create a secondary y-axis.

        Returns an opaque ``(ax, yaxis_name)`` handle for the ``*_secondary``
        primitives.
        """
        panel = _Panel.of(ax)
        secondary_y = panel.secondary_ref()
        panel.fig.update_layout(
            **{
                secondary_axis_key(secondary_y): dict(
                    overlaying=panel.ref("y"),
                    side="right",
                    anchor=panel.ref("x"),
                    showgrid=False,
                    showline=False,
                    zeroline=False,
                )
            }
        )

        return (ax, secondary_y)

    def line_secondary(
        self,
        secondary: Any,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        label: Optional[str] = None,
    ) -> Any:
        """Create a line plot on secondary y-axis."""
        ax, yaxis_name = secondary
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        dash = _DASH_MAP.get(linestyle, "solid")

        # For secondary axes, we need to set both xaxis and yaxis explicitly
        # and NOT use row/col which would override these references
        xaxis_ref = panel.ref("x")

        trace = go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color=color, width=linewidth, dash=dash),
            opacity=alpha,
            name=label or "",
            showlegend=label is not None,
            xaxis=xaxis_ref,
            yaxis=yaxis_name,
            hoverinfo="skip",
        )

        # Add trace directly without row/col to preserve axis references
        fig.add_trace(trace)
        return trace

    def fill_between_secondary(
        self,
        secondary: Any,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
    ) -> Any:
        """Fill area between two y-values on secondary y-axis."""
        ax, yaxis_name = secondary
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel

        if isinstance(y1, (int, float)):
            y1 = pd.Series([y1] * len(x))

        # For secondary axes, we need to set both xaxis and yaxis explicitly
        # and NOT use row/col which would override these references
        xaxis_ref = panel.ref("x")

        trace = go.Scatter(
            x=pd.concat([x, x[::-1]]),
            y=pd.concat([y2, y1[::-1]]),
            fill="toself",
            fillcolor=color,
            opacity=alpha,
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
            xaxis=xaxis_ref,
            yaxis=yaxis_name,
        )

        # Add trace directly without row/col to preserve axis references
        fig.add_trace(trace)
        return trace

    def set_secondary_ylim(
        self,
        secondary: Any,
        bottom: float,
        top: float,
    ) -> None:
        """Set secondary y-axis limits."""
        ax, yaxis_name = secondary
        _Panel.of(ax).fig.update_layout(
            **{secondary_axis_key(yaxis_name): dict(range=[bottom, top])}
        )

    def set_secondary_ylabel(
        self,
        secondary: Any,
        label: str,
        color: str = "black",
        fontsize: int = 10,
    ) -> None:
        """Set secondary y-axis label."""
        ax, yaxis_name = secondary
        label = convert_latex_to_unicode(label)
        _Panel.of(ax).fig.update_layout(
            **{
                secondary_axis_key(yaxis_name): dict(
                    title=dict(text=label, font=dict(size=fontsize, color=color)),
                    tickfont=dict(color=color),
                )
            }
        )

    def _add_legend_item(
        self,
        fig: go.Figure,
        row: int,
        name: str,
        color: str,
        symbol: str,
        size: int,
        legend_group: str,
        edgecolor: str = "black",
    ) -> None:
        """Add an invisible scatter trace for a legend entry."""
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    symbol=symbol,
                    size=size,
                    color=color,
                    line=dict(color=edgecolor, width=0.5),
                ),
                name=name,
                showlegend=True,
                legend=legend_group,
            ),
            row=row,
            col=1,
        )

    def add_panel_label(
        self,
        ax: Tuple[go.Figure, int],
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col
        fig.add_annotation(
            text=f"<b>{label}</b>",
            xref="x domain",
            yref="y domain",
            x=x_frac,
            y=y_frac,
            showarrow=False,
            font=dict(size=12),
            row=row,
            col=col,
        )

    def add_legend(
        self,
        ax: Tuple[go.Figure, int],
        entries: List[LegendEntry],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> None:
        """Render legend entries as an independently-positioned Plotly legend.

        Each call allocates a fresh legend key (legend, legend2, ...) so several
        legends coexist on one figure, positioned per panel row.
        """
        panel = _Panel.of(ax)
        fig, row = panel.fig, panel.row
        # to_plotly_json reports only legends that have been configured, so the
        # count of existing ones gives the next free key.
        existing = [k for k in fig.to_plotly_json()["layout"] if k.startswith("legend")]
        count = len(existing) + 1
        legend_key = "legend" if count == 1 else f"legend{count}"
        for entry in entries:
            symbol = (
                "square"
                if entry.marker == "patch"
                else _MARKER_SYMBOLS.get(entry.marker, "circle")
            )
            self._add_legend_item(
                fig,
                row,
                entry.label,
                entry.color,
                symbol,
                10,
                legend_key,
                entry.edgecolor or "black",
            )
        configure_legend(
            fig, row, legend_key, convert_latex_to_unicode(title or ""), loc
        )

    def hide_yaxis(self, ax: Tuple[go.Figure, int]) -> None:
        """Hide y-axis ticks, labels, line, and grid for gene track panels."""
        panel = _Panel.of(ax)
        fig, row, col, n_cols = panel
        fig.update_layout(
            **{
                panel.axis("yaxis"): dict(
                    showticklabels=False,
                    showline=False,
                    showgrid=False,
                    ticks="",
                )
            }
        )

    def format_xaxis_mb(self, ax: Tuple[go.Figure, int]) -> None:
        """Format x-axis to show megabase values."""
        panel = _Panel.of(ax)
        xaxis_name = panel.axis("xaxis")
        panel_range = x_range(panel, xaxis_name)
        if not panel_range:
            return
        tickvals, ticktext = mb_tick_positions(panel_range[0], panel_range[1])
        panel.fig.update_layout(
            **{
                xaxis_name: dict(
                    tickvals=tickvals,
                    ticktext=ticktext,
                    ticksuffix=" Mb",
                )
            }
        )

    def axvline(
        self,
        ax: Tuple[go.Figure, int],
        x: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a vertical line across the panel."""
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col
        dash = _DASH_MAP.get(linestyle, "dash")

        fig.add_vline(
            x=x,
            line_dash=dash,
            line_color=color,
            line_width=linewidth,
            opacity=alpha,
            row=row,
            col=col,
        )

    def hbar(
        self,
        ax: Tuple[go.Figure, int],
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
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        # Convert left to array if scalar
        if isinstance(left, (int, float)):
            left_arr = [left] * len(y)
        else:
            left_arr = list(left) if hasattr(left, "tolist") else left

        trace = go.Bar(
            y=y,
            x=width,
            orientation="h",
            base=left_arr,
            marker=dict(
                color=color,
                line=dict(color=edgecolor, width=linewidth),
            ),
            showlegend=False,
        )

        fig.add_trace(trace, row=row, col=col)
        return trace

    def errorbar_h(
        self,
        ax: Tuple[go.Figure, int],
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
        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        trace = go.Scatter(
            x=x,
            y=y,
            mode="markers",
            marker=dict(size=0),
            error_x=dict(
                type="data",
                symmetric=False,
                array=xerr_upper,
                arrayminus=xerr_lower,
                color=color,
                thickness=linewidth,
                width=capsize,
            ),
            showlegend=False,
        )

        fig.add_trace(trace, row=row, col=col)
        return trace

    def finalize_layout(
        self,
        fig: go.Figure,
        left: float = 0.08,
        right: float = 0.95,
        top: float = 0.95,
        bottom: float = 0.1,
        hspace: float = 0.08,
    ) -> None:
        """Adjust layout margins."""
        fig.update_layout(
            margin=dict(
                l=int(left * fig.layout.width) if fig.layout.width else 80,
                r=int((1 - right) * fig.layout.width) if fig.layout.width else 50,
                t=int((1 - top) * fig.layout.height) if fig.layout.height else 50,
                b=int(bottom * fig.layout.height) if fig.layout.height else 80,
            )
        )

    def add_region_highlight(
        self,
        fig: go.Figure,
        axes: List[Any],
        x_start: float,
        x_end: float,
        color: str = "yellow",
        alpha: float = 0.3,
    ) -> None:
        """Highlight an x-range across multiple plotly subplot rows."""
        for row in range(1, len(axes) + 1):
            fig.add_vrect(
                x0=x_start,
                x1=x_end,
                fillcolor=color,
                opacity=alpha,
                layer="below",
                line_width=0,
                row=row,
                col=1,
            )

    def highlight_heatmap_snp(
        self,
        ax: Tuple[go.Figure, int],
        fig: go.Figure,
        snp_idx: int,
        n_snps: int,
        color: str = "#FF0000",
        linewidth: float = 2,
    ) -> None:
        """Highlight a SNP's row and column in the heatmap."""
        for x, y in heatmap_highlight_cells(snp_idx, n_snps):
            fig.add_shape(
                type="rect",
                x0=x - 0.5,
                x1=x + 0.5,
                y0=y - 0.5,
                y1=y + 0.5,
                line=dict(color=color, width=linewidth),
                fillcolor="rgba(0,0,0,0)",
            )

    def add_heatmap(
        self,
        ax: Tuple[go.Figure, int],
        data: Any,
        x_coords: List[float],
        y_coords: List[float],
        cmap_colors: Optional[List[str]] = None,
        vmin: float = 0.0,
        vmax: float = 1.0,
    ) -> Any:
        """Render a heatmap of an already-shaped matrix."""
        import numpy as np

        panel = _Panel.of(ax)
        fig, row, col = panel.fig, panel.row, panel.col

        if cmap_colors is None:
            cmap_colors = ["#FFFFFF", "#FF0000"]

        # Plotly leaves a cell empty for None, not NaN.
        filled = np.ma.filled(np.ma.asarray(data).astype(float), np.nan)
        z = np.where(np.isnan(filled), None, filled)

        colorscale = [[0, cmap_colors[0]], [1, cmap_colors[1]]]

        fig.add_trace(
            go.Heatmap(
                z=z.tolist(),
                x=list(x_coords),
                y=list(y_coords),
                colorscale=colorscale,
                zmin=vmin,
                zmax=vmax,
                showscale=False,
            ),
            row=row,
            col=col,
        )
        # add_trace stores a copy, so hand back the figure's own trace: that is
        # the object add_colorbar has to mutate for the scale to appear.
        return fig.data[-1]

    def add_colorbar(
        self,
        ax: Tuple[go.Figure, int],
        mappable: Any,
        label: str = "R²",
        orientation: str = "vertical",
    ) -> Any:
        """Add colorbar legend for heatmap.

        Plotly draws the scale as part of the heatmap trace rather than as a
        separate artist, so this turns the trace's own scale on and titles it.
        ``add_heatmap`` leaves it off, which is what lets a caller skip this
        call to get a heatmap with no scale.
        """
        mappable.update(
            showscale=True,
            colorbar=dict(
                title=label,
                orientation="h" if orientation == "horizontal" else "v",
            ),
        )
        return mappable

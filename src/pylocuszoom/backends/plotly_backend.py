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
from .composition import LegendEntry, mb_tick_positions
from .hover import plotly_hovertemplate
from .plotly_layout import (
    _Panel,
    _SecondaryAxis,
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
_AXIS_STYLE = dict(
    showgrid=False,
    showline=True,
    linecolor="black",
    ticks="outside",
    zeroline=False,
)


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
    ) -> Tuple[go.Figure, List[_Panel]]:
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

        for row in range(1, n_panels + 1):
            panel = _Panel(fig, row)
            fig.update_layout(
                **{panel.axis("xaxis"): _AXIS_STYLE, panel.axis("yaxis"): _AXIS_STYLE}
            )

        return fig, [_Panel(fig, row) for row in range(1, n_panels + 1)]

    def create_figure_grid(
        self,
        n_rows: int,
        n_cols: int,
        width_ratios: Optional[List[float]] = None,
        height_ratios: Optional[List[float]] = None,
        figsize: Tuple[float, float] = (12.0, 8.0),
    ) -> Tuple[go.Figure, List[_Panel]]:
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

        for row in range(1, n_rows + 1):
            for col in range(1, n_cols + 1):
                panel = _Panel(fig, row, col, n_cols)
                fig.update_layout(
                    **{
                        panel.axis("xaxis"): _AXIS_STYLE,
                        panel.axis("yaxis"): _AXIS_STYLE,
                    }
                )

        return fig, [
            _Panel(fig, row, col, n_cols)
            for row in range(1, n_rows + 1)
            for col in range(1, n_cols + 1)
        ]

    def scatter(
        self,
        ax: _Panel,
        x: pd.Series,
        y: pd.Series,
        colors: Union[str, List[str], pd.Series],
        sizes: Union[float, List[float], pd.Series] = 60,
        marker: str = "o",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
        hover_data: Optional[pd.DataFrame] = None,
    ) -> Any:
        """Create a scatter plot on the given panel."""
        fig, row, col = ax.fig, ax.row, ax.col

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
            name="",
            showlegend=False,
        )

        fig.add_trace(trace, row=row, col=col)
        return trace

    def line(
        self,
        ax: Union[_Panel, _SecondaryAxis],
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        zorder: int = 1,
    ) -> Any:
        """Create a line plot on the given panel or secondary axis."""
        dash = _DASH_MAP.get(linestyle, "solid")
        # A secondary trace decorates the panel's data, so it skips the hover.
        hoverinfo = "skip" if isinstance(ax, _SecondaryAxis) else None

        trace = go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color=color, width=linewidth, dash=dash),
            opacity=alpha,
            name="",
            showlegend=False,
            xaxis=ax.xref,
            yaxis=ax.yref,
            hoverinfo=hoverinfo,
        )

        ax.fig.add_trace(trace)
        return trace

    def fill_between(
        self,
        ax: Union[_Panel, _SecondaryAxis],
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> Any:
        """Fill area between two y-values."""
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
            xaxis=ax.xref,
            yaxis=ax.yref,
        )

        ax.fig.add_trace(trace)
        return trace

    def axhline(
        self,
        ax: _Panel,
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the panel."""
        fig, row, col = ax.fig, ax.row, ax.col
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
        ax: _Panel,
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
        fig, row, col = ax.fig, ax.row, ax.col

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
        ax: _Panel,
        xy: Tuple[float, float],
        width: float,
        height: float,
        facecolor: Optional[str] = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a rectangle to the panel."""
        fig, row, col = ax.fig, ax.row, ax.col

        x0, y0 = xy
        x1, y1 = x0 + width, y0 + height

        fig.add_shape(
            type="rect",
            x0=x0,
            y0=y0,
            x1=x1,
            y1=y1,
            fillcolor=facecolor or "rgba(0,0,0,0)",
            line=dict(color=edgecolor, width=linewidth),
            row=row,
            col=col,
        )

    def add_polygon(
        self,
        ax: _Panel,
        points: List[List[float]],
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a polygon (e.g., triangle for strand arrows) to the panel."""
        fig, row, col = ax.fig, ax.row, ax.col

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

    def set_xlim(self, ax: _Panel, left: float, right: float) -> None:
        """Set x-axis limits."""
        ax.fig.update_layout(**{ax.axis("xaxis"): dict(range=[left, right])})

    def set_ylim(self, ax: _Panel, bottom: float, top: float) -> None:
        """Set y-axis limits."""
        ax.fig.update_layout(**{ax.axis("yaxis"): dict(range=[bottom, top])})

    def set_xlabel(self, ax: _Panel, label: str, fontsize: int = 12) -> None:
        """Set x-axis label."""
        label = convert_latex_to_unicode(label)
        ax.fig.update_layout(
            **{ax.axis("xaxis"): dict(title=dict(text=label, font=dict(size=fontsize)))}
        )

    def set_ylabel(self, ax: _Panel, label: str, fontsize: int = 12) -> None:
        """Set y-axis label."""
        label = convert_latex_to_unicode(label)
        ax.fig.update_layout(
            **{ax.axis("yaxis"): dict(title=dict(text=label, font=dict(size=fontsize)))}
        )

    def set_yticks(
        self,
        ax: _Panel,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
    ) -> None:
        """Set y-axis tick positions and labels."""
        ax.fig.update_layout(
            **{
                ax.axis("yaxis"): dict(
                    tickmode="array",
                    tickvals=positions,
                    ticktext=labels,
                    tickfont=dict(size=fontsize),
                )
            }
        )

    def set_xticks(
        self,
        ax: _Panel,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
        rotation: int = 0,
        ha: str = "center",
    ) -> None:
        """Set x-axis tick positions and labels."""
        ax.fig.update_layout(
            **{
                ax.axis("xaxis"): dict(
                    tickmode="array",
                    tickvals=positions,
                    ticktext=labels,
                    tickfont=dict(size=fontsize),
                    tickangle=-rotation if rotation else 0,
                )
            }
        )

    def set_title(self, ax: _Panel, title: str, fontsize: int = 14) -> None:
        """Set subplot title using annotation.

        For grid layouts, this adds an annotation above the subplot.
        For single-column layouts, sets the global figure title for the first panel.
        """
        if ax.n_cols == 1 and ax.row == 1:
            # Single-column layout: use global figure title
            ax.fig.update_layout(title=dict(text=title, font=dict(size=fontsize)))
        else:
            # Grid layout: add annotation above the subplot
            # Use subplot's axis domain for positioning
            xref = f"{ax.ref('x')} domain"
            yref = f"{ax.ref('y')} domain"

            ax.fig.add_annotation(
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

    def create_twin_axis(self, ax: _Panel) -> _SecondaryAxis:
        """Create a secondary y-axis and return its handle."""
        secondary_y = ax.secondary_ref()
        ax.fig.update_layout(
            **{
                secondary_axis_key(secondary_y): dict(
                    overlaying=ax.ref("y"),
                    side="right",
                    anchor=ax.ref("x"),
                    showgrid=False,
                    showline=False,
                    zeroline=False,
                )
            }
        )

        return _SecondaryAxis(ax, secondary_y)

    def set_secondary_ylim(
        self,
        secondary: _SecondaryAxis,
        bottom: float,
        top: float,
    ) -> None:
        """Set secondary y-axis limits."""
        secondary.fig.update_layout(
            **{secondary_axis_key(secondary.yref): dict(range=[bottom, top])}
        )

    def set_secondary_ylabel(
        self,
        secondary: _SecondaryAxis,
        label: str,
        color: str = "black",
        fontsize: int = 10,
    ) -> None:
        """Set secondary y-axis label."""
        label = convert_latex_to_unicode(label)
        secondary.fig.update_layout(
            **{
                secondary_axis_key(secondary.yref): dict(
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
        ax: _Panel,
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        fig, row, col = ax.fig, ax.row, ax.col
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
        ax: _Panel,
        entries: List[LegendEntry],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> None:
        """Render legend entries as an independently-positioned Plotly legend.

        Each call allocates a fresh legend key (legend, legend2, ...) so several
        legends coexist on one figure, positioned per panel row.
        """
        fig, row = ax.fig, ax.row
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

    def hide_yaxis(self, ax: _Panel) -> None:
        """Hide y-axis ticks, labels, line, and grid for gene track panels."""
        ax.fig.update_layout(
            **{
                ax.axis("yaxis"): dict(
                    showticklabels=False,
                    showline=False,
                    showgrid=False,
                    ticks="",
                )
            }
        )

    def format_xaxis_mb(self, ax: _Panel) -> None:
        """Format x-axis to show megabase values."""
        xaxis_name = ax.axis("xaxis")
        panel_range = x_range(ax, xaxis_name)
        if not panel_range:
            return
        tickvals, ticktext = mb_tick_positions(panel_range[0], panel_range[1])
        ax.fig.update_layout(
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
        ax: _Panel,
        x: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a vertical line across the panel."""
        fig, row, col = ax.fig, ax.row, ax.col
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

    def errorbar_h(
        self,
        ax: _Panel,
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
        fig, row, col = ax.fig, ax.row, ax.col

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

    def add_heatmap(
        self,
        ax: _Panel,
        data: Any,
        x_coords: List[float],
        y_coords: List[float],
        cmap_colors: List[str],
        vmin: float = 0.0,
        vmax: float = 1.0,
    ) -> Any:
        """Render a heatmap of an already-shaped matrix."""
        import numpy as np

        fig, row, col = ax.fig, ax.row, ax.col

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
        ax: _Panel,
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

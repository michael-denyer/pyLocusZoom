"""Plotly backend for pyLocusZoom.

Interactive backend with hover tooltips and zoom/pan capabilities.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from . import convert_latex_to_unicode, register_backend
from .composition import LegendEntry, heatmap_highlight_cells
from .hover import plotly_hovertemplate

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
# Matplotlib legend `loc` vocabulary mapped to Plotly (xanchor, yanchor).
# "best" has no Plotly equivalent, so it takes the upper-right default.
_LEGEND_ANCHORS = {
    "best": ("right", "top"),
    "upper right": ("right", "top"),
    "upper left": ("left", "top"),
    "upper center": ("center", "top"),
    "lower right": ("right", "bottom"),
    "lower left": ("left", "bottom"),
    "lower center": ("center", "bottom"),
    "center right": ("right", "middle"),
    "center left": ("left", "middle"),
    "right": ("right", "middle"),
    "center": ("center", "middle"),
}
_LEGEND_X = {"left": 0.01, "center": 0.5, "right": 0.99}


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
        """Create a figure with multiple panels.

        Args:
            n_panels: Number of vertical panels.
            height_ratios: Relative heights for each panel.
            figsize: Figure size as (width, height) in inches.
            sharex: Whether panels share the x-axis.

        Returns:
            Tuple of (figure, list of row indices for each panel).
        """
        # Convert inches to pixels (assuming 100 dpi for web)
        width_px = int(figsize[0] * 100)
        height_px = int(figsize[1] * 100)

        # Normalize height ratios
        total = sum(height_ratios)
        row_heights = [h / total for h in height_ratios]

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
            xaxis = self._axis_name("xaxis", row)
            yaxis = self._axis_name("yaxis", row)
            fig.update_layout(**{xaxis: axis_style, yaxis: axis_style})

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
        """Create a figure with a grid of subplots.

        Args:
            n_rows: Number of rows.
            n_cols: Number of columns.
            width_ratios: Relative widths for columns.
            height_ratios: Relative heights for rows.
            figsize: Figure size as (width, height).

        Returns:
            Tuple of (figure, flattened list of (fig, row, col) tuples).
        """
        width_px = int(figsize[0] * 100)
        height_px = int(figsize[1] * 100)

        # Normalize ratios
        if width_ratios is not None:
            total = sum(width_ratios)
            column_widths = [w / total for w in width_ratios]
        else:
            column_widths = None

        if height_ratios is not None:
            total = sum(height_ratios)
            row_heights_norm = [h / total for h in height_ratios]
        else:
            row_heights_norm = None

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
                subplot_idx = (row - 1) * n_cols + col
                xaxis = f"xaxis{subplot_idx}" if subplot_idx > 1 else "xaxis"
                yaxis = f"yaxis{subplot_idx}" if subplot_idx > 1 else "yaxis"
                fig.update_layout(**{xaxis: axis_style, yaxis: axis_style})

        # Return flattened list of (fig, row, col, n_cols) tuples
        panel_refs = []
        for row in range(1, n_rows + 1):
            for col in range(1, n_cols + 1):
                panel_refs.append((fig, row, col, n_cols))
        return fig, panel_refs

    def _extract_row_col(self, ax: Any) -> Tuple[go.Figure, int, int, int]:
        """Extract figure, row, col, and n_cols from ax tuple.

        Handles (fig, row), (fig, row, col), and (fig, row, col, n_cols) formats.
        The 3-tuple format assumes single-column layout (n_cols=1).

        Returns:
            Tuple of (figure, row, col, n_cols).

        Raises:
            ValueError: If ax tuple has unexpected length.
        """
        if len(ax) == 2:
            fig, row = ax
            return fig, row, 1, 1
        elif len(ax) == 3:
            fig, row, col = ax
            return fig, row, col, 1
        elif len(ax) == 4:
            fig, row, col, n_cols = ax
            return fig, row, col, n_cols
        else:
            raise ValueError(f"Expected ax tuple of length 2, 3, or 4, got {len(ax)}")

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
        fig, row, col, _ = self._extract_row_col(ax)

        # Convert matplotlib marker to plotly symbol
        symbol = _MARKER_SYMBOLS.get(marker, "circle")

        # Convert size (matplotlib uses area, plotly uses diameter)
        if isinstance(sizes, (int, float)):
            size = max(6, sizes**0.5)  # Approximate conversion
        else:
            size = [max(6, s**0.5) for s in sizes]

        # Build hover template
        if hover_data is not None:
            customdata = hover_data.values
            hovertemplate = plotly_hovertemplate(hover_data)
        else:
            customdata = None
            hovertemplate = "x: %{x}<br>y: %{y:.2f}<extra></extra>"

        # Handle color - could be single color or array
        if isinstance(colors, str):
            marker_color = colors
        else:
            marker_color = list(colors) if hasattr(colors, "tolist") else colors

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
        fig, row, col, _ = self._extract_row_col(ax)
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
        fig, row, col, _ = self._extract_row_col(ax)

        # Convert y1 to series if scalar
        if isinstance(y1, (int, float)):
            y1 = pd.Series([y1] * len(x))

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
        fig, row, col, _ = self._extract_row_col(ax)
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
        fig, row, col, _ = self._extract_row_col(ax)

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
        fig, row, col, _ = self._extract_row_col(ax)

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
        fig, row, col, _ = self._extract_row_col(ax)

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
        fig, row, col, n_cols = self._extract_row_col(ax)
        fig.update_layout(
            **{self._axis_name("xaxis", row, col, n_cols): dict(range=[left, right])}
        )

    def set_ylim(self, ax: Tuple[go.Figure, int], bottom: float, top: float) -> None:
        """Set y-axis limits."""
        fig, row, col, n_cols = self._extract_row_col(ax)
        fig.update_layout(
            **{self._axis_name("yaxis", row, col, n_cols): dict(range=[bottom, top])}
        )

    def set_xlabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set x-axis label."""
        fig, row, col, n_cols = self._extract_row_col(ax)
        label = convert_latex_to_unicode(label)
        fig.update_layout(
            **{
                self._axis_name("xaxis", row, col, n_cols): dict(
                    title=dict(text=label, font=dict(size=fontsize))
                )
            }
        )

    def set_ylabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set y-axis label."""
        fig, row, col, n_cols = self._extract_row_col(ax)
        label = convert_latex_to_unicode(label)
        fig.update_layout(
            **{
                self._axis_name("yaxis", row, col, n_cols): dict(
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
        fig, row, col, n_cols = self._extract_row_col(ax)
        fig.update_layout(
            **{
                self._axis_name("yaxis", row, col, n_cols): dict(
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
        fig, row, col, n_cols = self._extract_row_col(ax)
        fig.update_layout(
            **{
                self._axis_name("xaxis", row, col, n_cols): dict(
                    tickmode="array",
                    tickvals=positions,
                    ticktext=labels,
                    tickfont=dict(size=fontsize),
                    tickangle=-rotation if rotation else 0,
                )
            }
        )

    def _axis_name(self, axis: str, row: int, col: int = 1, n_cols: int = 1) -> str:
        """Get Plotly axis name for a given row and column.

        Plotly names axes using a linear subplot index:
        - subplot (1,1) uses 'xaxis', 'yaxis'
        - subplot (1,2) uses 'xaxis2', 'yaxis2'
        - subplot (2,1) uses 'xaxis3', 'yaxis3' (for 2-column grid)

        Args:
            axis: Base axis name ('xaxis' or 'yaxis').
            row: Row number (1-indexed).
            col: Column number (1-indexed).
            n_cols: Total number of columns in the grid.

        Returns:
            Plotly axis name string.
        """
        subplot_idx = (row - 1) * n_cols + col
        return f"{axis}{subplot_idx}" if subplot_idx > 1 else axis

    def set_title(
        self, ax: Tuple[go.Figure, int], title: str, fontsize: int = 14
    ) -> None:
        """Set subplot title using annotation.

        For grid layouts, this adds an annotation above the subplot.
        For single-column layouts, sets the global figure title for the first panel.
        """
        fig, row, col, n_cols = self._extract_row_col(ax)

        if n_cols == 1 and row == 1:
            # Single-column layout: use global figure title
            fig.update_layout(title=dict(text=title, font=dict(size=fontsize)))
        else:
            # Grid layout: add annotation above the subplot
            # Use subplot's axis domain for positioning
            # Plotly uses "x" or "x2", "x3", etc. (not "xaxis")
            subplot_idx = (row - 1) * n_cols + col
            xref = f"x{subplot_idx} domain" if subplot_idx > 1 else "x domain"
            yref = f"y{subplot_idx} domain" if subplot_idx > 1 else "y domain"

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

        For Plotly subplots, we need unique axis names that don't conflict
        with the subplot axes. We use a high number suffix to avoid conflicts.
        """
        fig, row, col, n_cols = self._extract_row_col(ax)

        # Calculate subplot index for proper axis naming
        subplot_idx = (row - 1) * n_cols + col

        # Use a high offset to avoid collision with primary subplot axes.
        # Plotly names primary axes yaxis, yaxis2, ..., yaxisN for N subplots.
        # Offset of 100 supports up to 99 subplot rows without collision.
        secondary_suffix = 100 + subplot_idx - 1
        secondary_y = f"y{secondary_suffix}"
        yaxis_name = f"yaxis{secondary_suffix}"

        # Get the primary y-axis name for this subplot
        primary_y = f"y{subplot_idx}" if subplot_idx > 1 else "y"
        # Get the x-axis name for this subplot
        xaxis_ref = f"x{subplot_idx}" if subplot_idx > 1 else "x"

        # Configure secondary y-axis to overlay the primary axis of this row
        fig.update_layout(
            **{
                yaxis_name: dict(
                    overlaying=primary_y,
                    side="right",
                    anchor=xaxis_ref,
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
        fig, row, col, n_cols = self._extract_row_col(ax)
        dash = _DASH_MAP.get(linestyle, "solid")

        # For secondary axes, we need to set both xaxis and yaxis explicitly
        # and NOT use row/col which would override these references
        subplot_idx = (row - 1) * n_cols + col
        xaxis_ref = f"x{subplot_idx}" if subplot_idx > 1 else "x"

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
        fig, row, col, n_cols = self._extract_row_col(ax)

        if isinstance(y1, (int, float)):
            y1 = pd.Series([y1] * len(x))

        # For secondary axes, we need to set both xaxis and yaxis explicitly
        # and NOT use row/col which would override these references
        subplot_idx = (row - 1) * n_cols + col
        xaxis_ref = f"x{subplot_idx}" if subplot_idx > 1 else "x"

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
        fig, row, col, _ = self._extract_row_col(ax)
        yaxis_key = (
            "yaxis" + yaxis_name[1:] if yaxis_name.startswith("y") else yaxis_name
        )
        fig.update_layout(**{yaxis_key: dict(range=[bottom, top])})

    def set_secondary_ylabel(
        self,
        secondary: Any,
        label: str,
        color: str = "black",
        fontsize: int = 10,
    ) -> None:
        """Set secondary y-axis label."""
        ax, yaxis_name = secondary
        fig, row, col, _ = self._extract_row_col(ax)
        label = convert_latex_to_unicode(label)
        yaxis_key = (
            "yaxis" + yaxis_name[1:] if yaxis_name.startswith("y") else yaxis_name
        )
        fig.update_layout(
            **{
                yaxis_key: dict(
                    title=dict(text=label, font=dict(size=fontsize, color=color)),
                    tickfont=dict(color=color),
                )
            }
        )

    def _get_panel_y(self, fig: go.Figure, row: int, vertical: str) -> float:
        """Get a y-coordinate (in paper coords) within a subplot row's domain.

        Plotly subplots have y-axis domains that define their vertical position.
        ``vertical`` selects the top, bottom, or midpoint of that domain.
        """
        yaxis = getattr(fig.layout, self._axis_name("yaxis", row), None)
        domain = yaxis.domain if yaxis and yaxis.domain else (0.01, 0.99)
        if vertical == "bottom":
            return domain[0]
        if vertical == "middle":
            return (domain[0] + domain[1]) / 2
        return domain[1]

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

    def _configure_legend(
        self, fig: go.Figure, row: int, legend_key: str, title: str, loc: str
    ) -> None:
        """Configure legend position and styling."""
        horizontal, vertical = _LEGEND_ANCHORS.get(loc, _LEGEND_ANCHORS["upper right"])
        fig.update_layout(
            **{
                legend_key: dict(
                    title=dict(text=convert_latex_to_unicode(title)),
                    x=_LEGEND_X[horizontal],
                    y=self._get_panel_y(fig, row, vertical),
                    xanchor=horizontal,
                    yanchor=vertical,
                    bgcolor="rgba(255,255,255,0.9)",
                    bordercolor="black",
                    borderwidth=1,
                )
            }
        )

    def add_panel_label(
        self,
        ax: Tuple[go.Figure, int],
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        fig, row, col, _ = self._extract_row_col(ax)
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
        fig, row, col, _ = self._extract_row_col(ax)
        count = getattr(fig, "_legend_count", 0) + 1
        fig._legend_count = count
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
        self._configure_legend(fig, row, legend_key, title or "", loc)

    def hide_spines(self, ax: Tuple[go.Figure, int], spines: List[str]) -> None:
        """Hide specified axis spines (lines).

        Plotly doesn't have spines, but we can hide axis lines.
        """
        # Plotly's template "plotly_white" already hides top/right lines
        # No action needed - method exists for API compatibility
        pass

    def hide_yaxis(self, ax: Tuple[go.Figure, int]) -> None:
        """Hide y-axis ticks, labels, line, and grid for gene track panels."""
        fig, row, col, n_cols = self._extract_row_col(ax)
        fig.update_layout(
            **{
                self._axis_name("yaxis", row, col, n_cols): dict(
                    showticklabels=False,
                    showline=False,
                    showgrid=False,
                    ticks="",
                )
            }
        )

    def format_xaxis_mb(self, ax: Tuple[go.Figure, int]) -> None:
        """Format x-axis to show megabase values.

        Stores the subplot info for later tick formatting in finalize_layout.
        """
        fig, row, col, n_cols = self._extract_row_col(ax)
        # Store that this axis needs Mb formatting
        if not hasattr(fig, "_mb_format_rows"):
            fig._mb_format_rows = []
        # Store (row, col, n_cols) tuple for proper axis naming later
        fig._mb_format_rows.append((row, col, n_cols))

    def save(
        self,
        fig: go.Figure,
        path: str,
        dpi: int = 150,
        bbox_inches: str = "tight",
    ) -> None:
        """Save figure to file.

        Supports .html for interactive and .png/.pdf for static.
        """
        if path.endswith(".html"):
            fig.write_html(path)
        else:
            # Static export requires kaleido
            scale = dpi / 100
            fig.write_image(path, scale=scale)

    def show(self, fig: go.Figure) -> None:
        """Display the figure."""
        fig.show()

    def close(self, fig: go.Figure) -> None:
        """Close the figure (no-op for plotly)."""
        pass

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
        fig, row, col, _ = self._extract_row_col(ax)
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
        fig, row, col, _ = self._extract_row_col(ax)

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
        fig, row, col, _ = self._extract_row_col(ax)

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
        """Adjust layout margins and apply Mb tick formatting.

        Args:
            fig: Figure object.
            left, right, top, bottom: Margins as fractions.
            hspace: Ignored for plotly (use vertical_spacing in make_subplots).
        """
        fig.update_layout(
            margin=dict(
                l=int(left * fig.layout.width) if fig.layout.width else 80,
                r=int((1 - right) * fig.layout.width) if fig.layout.width else 50,
                t=int((1 - top) * fig.layout.height) if fig.layout.height else 50,
                b=int(bottom * fig.layout.height) if fig.layout.height else 80,
            )
        )

        # Apply Mb tick formatting to marked axes
        if hasattr(fig, "_mb_format_rows"):
            import numpy as np

            for item in fig._mb_format_rows:
                # Handle both old format (row) and new format (row, col, n_cols)
                if isinstance(item, tuple):
                    row, col, n_cols = item
                else:
                    row, col, n_cols = item, 1, 1
                xaxis_name = self._axis_name("xaxis", row, col, n_cols)
                xaxis = getattr(fig.layout, xaxis_name, None)

                # Get x-range from the axis or compute from data
                x_range = None
                if xaxis and xaxis.range:
                    x_range = xaxis.range
                else:
                    # Compute from trace data (filter out None values from legend traces)
                    x_vals = []
                    for trace in fig.data:
                        if hasattr(trace, "x") and trace.x is not None:
                            x_vals.extend([v for v in trace.x if v is not None])
                    if x_vals:
                        x_range = [min(x_vals), max(x_vals)]

                if x_range:
                    x_min_mb, x_max_mb = x_range[0] / 1e6, x_range[1] / 1e6
                    span_mb = x_max_mb - x_min_mb

                    # Choose tick spacing based on range
                    if span_mb <= 0.5:
                        tick_step = 0.1
                    elif span_mb <= 2:
                        tick_step = 0.25
                    elif span_mb <= 5:
                        tick_step = 0.5
                    elif span_mb <= 20:
                        tick_step = 2
                    else:
                        tick_step = 5

                    # Generate ticks
                    first_tick = np.ceil(x_min_mb / tick_step) * tick_step
                    tickvals_mb = np.arange(
                        first_tick, x_max_mb + tick_step / 2, tick_step
                    )
                    tickvals_bp = [v * 1e6 for v in tickvals_mb]
                    ticktext = [f"{v:.2f}" for v in tickvals_mb]

                    fig.update_layout(
                        **{
                            xaxis_name: dict(
                                tickvals=tickvals_bp,
                                ticktext=ticktext,
                                ticksuffix=" Mb",
                            )
                        }
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
        """Highlight a SNP's row and column in the heatmap.

        Args:
            ax: Tuple of (figure, row_number) (unused, shapes added to fig).
            fig: Plotly figure for adding shapes.
            snp_idx: Index of SNP to highlight.
            n_snps: Total number of SNPs in matrix.
            color: Highlight color.
            linewidth: Line width for highlight rectangles.
        """
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
        mask_upper: bool = True,
    ) -> Any:
        """Render heatmap with optional triangular masking.

        Args:
            ax: Tuple of (figure, row_number).
            data: 2D numpy array of values (NaN for missing).
            x_coords: X coordinates for cells.
            y_coords: Y coordinates for cells.
            cmap_colors: Color gradient endpoints [start, end].
            vmin: Minimum value for color scale.
            vmax: Maximum value for color scale.
            mask_upper: If True, mask upper triangle.

        Returns:
            Heatmap trace.
        """
        import numpy as np

        fig, row, col, _ = self._extract_row_col(ax)

        if cmap_colors is None:
            cmap_colors = ["#FFFFFF", "#FF0000"]

        # Mask upper triangle by setting to NaN
        plot_data = data.copy()
        if mask_upper:
            for i in range(data.shape[0]):
                for j in range(i + 1, data.shape[1]):
                    plot_data[i, j] = np.nan

        # Replace NaN with None for Plotly
        z = np.where(np.isnan(plot_data), None, plot_data)

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

        Args:
            ax: Tuple of (figure, row_number) (unused, the trace carries it).
            mappable: Heatmap trace returned by add_heatmap.
            label: Colorbar label.
            orientation: "vertical" or "horizontal".

        Returns:
            The heatmap trace whose scale was enabled.
        """
        mappable.update(
            showscale=True,
            colorbar=dict(
                title=label,
                orientation="h" if orientation == "horizontal" else "v",
            ),
        )
        return mappable

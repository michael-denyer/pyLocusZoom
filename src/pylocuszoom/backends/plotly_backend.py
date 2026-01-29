"""Plotly backend for pyLocusZoom.

Interactive backend with hover tooltips and zoom/pan capabilities.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from . import convert_latex_to_unicode, register_backend

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
    def supports_snp_labels(self) -> bool:
        """Plotly uses hover tooltips instead of labels."""
        return False

    @property
    def supports_hover(self) -> bool:
        """Plotly supports hover tooltips."""
        return True

    @property
    def supports_secondary_axis(self) -> bool:
        """Plotly supports secondary y-axis."""
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

        # Store grid dimensions for axis naming
        fig._n_cols = n_cols
        fig._n_rows = n_rows

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

        # Return flattened list of (fig, row, col) tuples
        panel_refs = []
        for row in range(1, n_rows + 1):
            for col in range(1, n_cols + 1):
                panel_refs.append((fig, row, col))
        return fig, panel_refs

    def _extract_row_col(self, ax: Any) -> Tuple[go.Figure, int, int, int]:
        """Extract figure, row, col, and n_cols from ax tuple.

        Handles both (fig, row) for create_figure and (fig, row, col) for create_figure_grid.

        Returns:
            Tuple of (figure, row, col, n_cols).
        """
        if len(ax) == 2:
            fig, row = ax
            col = 1
        else:
            fig, row, col = ax
        # Get n_cols from figure if stored, default to 1 for single-column layouts
        n_cols = getattr(fig, "_n_cols", 1)
        return fig, row, col, n_cols

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
            hover_cols = hover_data.columns.tolist()
            hovertemplate = "<b>%{customdata[0]}</b><br>"
            for i, col_name in enumerate(hover_cols[1:], 1):
                col_lower = col_name.lower()
                if col_lower in ("p-value", "pval", "p_value"):
                    hovertemplate += f"{col_name}: %{{customdata[{i}]:.2e}}<br>"
                elif any(x in col_lower for x in ("r2", "r²", "ld")):
                    hovertemplate += f"{col_name}: %{{customdata[{i}]:.3f}}<br>"
                elif "pos" in col_lower:
                    hovertemplate += f"{col_name}: %{{customdata[{i}]:,.0f}}<br>"
                else:
                    hovertemplate += f"{col_name}: %{{customdata[{i}]}}<br>"
            hovertemplate += "<extra></extra>"
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
        label = self._convert_label(label)
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
        label = self._convert_label(label)
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

    def _get_legend_position(self, loc: str) -> dict:
        """Map matplotlib-style legend location to Plotly position dict."""
        loc_map = {
            "upper left": dict(x=0.01, y=0.99, xanchor="left", yanchor="top"),
            "upper right": dict(x=0.99, y=0.99, xanchor="right", yanchor="top"),
            "lower left": dict(x=0.01, y=0.01, xanchor="left", yanchor="bottom"),
            "lower right": dict(x=0.99, y=0.01, xanchor="right", yanchor="bottom"),
        }
        return loc_map.get(loc, loc_map["upper left"])

    def _convert_label(self, label: str) -> str:
        """Convert LaTeX-style labels to Unicode for Plotly display."""
        return convert_latex_to_unicode(label)

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

    def create_twin_axis(self, ax: Tuple[go.Figure, int]) -> Tuple[go.Figure, int, str]:
        """Create a secondary y-axis.

        Returns tuple of (figure, row, secondary_yaxis_name).

        For Plotly subplots, we need unique axis names that don't conflict
        with the subplot axes. We use a high number suffix to avoid conflicts.
        """
        fig, row, col, n_cols = self._extract_row_col(ax)

        # Calculate subplot index for proper axis naming
        subplot_idx = (row - 1) * n_cols + col

        # Use a unique suffix that won't conflict with subplot axis numbering
        # yaxis10, yaxis11, etc. are unlikely to conflict with typical subplot counts
        secondary_suffix = 10 + subplot_idx - 1
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

        return (fig, row, secondary_y)

    def line_secondary(
        self,
        ax: Tuple[go.Figure, int],
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        label: Optional[str] = None,
        yaxis_name: str = "y2",
    ) -> Any:
        """Create a line plot on secondary y-axis."""
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
        ax: Tuple[go.Figure, int],
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        yaxis_name: str = "y2",
    ) -> Any:
        """Fill area between two y-values on secondary y-axis."""
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
        ax: Tuple[go.Figure, int],
        bottom: float,
        top: float,
        yaxis_name: str = "y2",
    ) -> None:
        """Set secondary y-axis limits."""
        fig, row, col, _ = self._extract_row_col(ax)
        yaxis_key = (
            "yaxis" + yaxis_name[1:] if yaxis_name.startswith("y") else yaxis_name
        )
        fig.update_layout(**{yaxis_key: dict(range=[bottom, top])})

    def set_secondary_ylabel(
        self,
        ax: Tuple[go.Figure, int],
        label: str,
        color: str = "black",
        fontsize: int = 10,
        yaxis_name: str = "y2",
    ) -> None:
        """Set secondary y-axis label."""
        fig, row, col, _ = self._extract_row_col(ax)
        label = self._convert_label(label)
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

    def _get_panel_y_top(self, fig: go.Figure, row: int) -> float:
        """Get the top y-coordinate (in paper coords) for a subplot row.

        Plotly subplots have y-axis domains that define their vertical position.
        This returns the top of the domain for positioning legends.
        """
        yaxis = getattr(fig.layout, self._axis_name("yaxis", row), None)
        if yaxis and yaxis.domain:
            return yaxis.domain[1]
        return 0.99

    def _add_legend_item(
        self,
        fig: go.Figure,
        row: int,
        name: str,
        color: str,
        symbol: str,
        size: int,
        legend_group: str,
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
                    line=dict(color="black", width=0.5),
                ),
                name=name,
                showlegend=True,
                legend=legend_group,
            ),
            row=row,
            col=1,
        )

    def _configure_legend(
        self, fig: go.Figure, row: int, legend_key: str, title: str
    ) -> None:
        """Configure legend position and styling."""
        y_pos = self._get_panel_y_top(fig, row)
        fig.update_layout(
            **{
                legend_key: dict(
                    title=dict(text=title),
                    x=0.99,
                    y=y_pos,
                    xanchor="right",
                    yanchor="top",
                    bgcolor="rgba(255,255,255,0.9)",
                    bordercolor="black",
                    borderwidth=1,
                )
            }
        )

    def add_snp_labels(
        self,
        ax: Tuple[go.Figure, int],
        df: pd.DataFrame,
        pos_col: str,
        neglog10p_col: str,
        rs_col: str,
        label_top_n: int,
        genes_df: Optional[pd.DataFrame],
        chrom: int,
    ) -> None:
        """No-op: Plotly uses hover tooltips instead of text labels."""
        pass

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

    def add_ld_legend(
        self,
        ax: Tuple[go.Figure, int],
        ld_bins: List[Tuple[float, str, str]],
        lead_snp_color: str,
    ) -> None:
        """Add LD color legend using invisible scatter traces.

        Uses Plotly's separate legend feature (legend="legend") so LD legend
        can be positioned independently from eQTL and fine-mapping legends.
        """
        fig, row, col, _ = self._extract_row_col(ax)

        self._add_legend_item(
            fig, row, "Lead SNP", lead_snp_color, "diamond", 12, "legend"
        )
        for _, label, color in ld_bins:
            self._add_legend_item(fig, row, label, color, "square", 10, "legend")

        self._configure_legend(fig, row, "legend", "r²")

    def add_legend(
        self,
        ax: Tuple[go.Figure, int],
        handles: List[Any],
        labels: List[str],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> Any:
        """Add a legend to the figure.

        Note: Plotly handles legends automatically from trace names.
        This method updates legend positioning.
        """
        fig, _ = ax
        legend_pos = self._get_legend_position(loc)
        fig.update_layout(
            legend=dict(
                **legend_pos,
                title=dict(text=title) if title else None,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        )

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

    def add_eqtl_legend(
        self,
        ax: Tuple[go.Figure, int],
        eqtl_positive_bins: List[Tuple[float, float, str, str]],
        eqtl_negative_bins: List[Tuple[float, float, str, str]],
    ) -> None:
        """Add eQTL effect size legend using invisible scatter traces.

        Uses Plotly's separate legend feature (legend="legend2") so eQTL legend
        is positioned independently below the LD legend.
        """
        fig, row, col, _ = self._extract_row_col(ax)

        for _, _, label, color in eqtl_positive_bins:
            self._add_legend_item(fig, row, label, color, "triangle-up", 10, "legend2")
        for _, _, label, color in eqtl_negative_bins:
            self._add_legend_item(
                fig, row, label, color, "triangle-down", 10, "legend2"
            )

        self._configure_legend(fig, row, "legend2", "eQTL effect")

    def add_finemapping_legend(
        self,
        ax: Tuple[go.Figure, int],
        credible_sets: List[int],
        get_color_func: Any,
    ) -> None:
        """Add fine-mapping credible set legend using invisible scatter traces.

        Uses Plotly's separate legend feature (legend="legend2") so fine-mapping
        legend is positioned independently below the LD legend.
        """
        if not credible_sets:
            return

        fig, row, col, _ = self._extract_row_col(ax)

        for cs_id in credible_sets:
            self._add_legend_item(
                fig, row, f"CS{cs_id}", get_color_func(cs_id), "circle", 10, "legend2"
            )

        self._configure_legend(fig, row, "legend2", "Credible sets")

    def add_simple_legend(
        self,
        ax: Tuple[go.Figure, int],
        label: str,
        loc: str = "upper right",
    ) -> None:
        """Add simple legend positioning.

        Plotly handles legends automatically from trace names.
        This just positions the legend.
        """
        fig, _ = ax
        legend_pos = self._get_legend_position(loc)
        fig.update_layout(
            legend=dict(
                **legend_pos,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
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

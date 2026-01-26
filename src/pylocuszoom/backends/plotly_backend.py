"""Plotly backend for pyLocusZoom.

Interactive backend with hover tooltips and zoom/pan capabilities.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


class PlotlyBackend:
    """Plotly backend for interactive plot generation.

    Produces interactive HTML plots with hover tooltips showing:
    - SNP RS ID
    - P-value
    - R² with lead SNP
    - Nearest gene
    """

    def __init__(self) -> None:
        """Initialize the plotly backend."""
        self._marker_symbols = {
            "o": "circle",
            "D": "diamond",
            "s": "square",
            "^": "triangle-up",
            "v": "triangle-down",
        }

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

        # Return row indices (1-based for plotly)
        panel_refs = list(range(1, n_panels + 1))
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

        For plotly, ax is a tuple of (figure, row_number).
        """
        fig, row = ax

        # Convert matplotlib marker to plotly symbol
        symbol = self._marker_symbols.get(marker, "circle")

        # Convert size (matplotlib uses area, plotly uses diameter)
        if isinstance(sizes, (int, float)):
            size = max(6, sizes ** 0.5)  # Approximate conversion
        else:
            size = [max(6, s ** 0.5) for s in sizes]

        # Build hover template
        if hover_data is not None:
            customdata = hover_data.values
            hover_cols = hover_data.columns.tolist()
            hovertemplate = "<b>%{customdata[0]}</b><br>"
            for i, col in enumerate(hover_cols[1:], 1):
                if "p" in col.lower():
                    hovertemplate += f"{col}: %{{customdata[{i}]:.2e}}<br>"
                elif "r2" in col.lower() or "ld" in col.lower():
                    hovertemplate += f"{col}: %{{customdata[{i}]:.3f}}<br>"
                else:
                    hovertemplate += f"{col}: %{{customdata[{i}]}}<br>"
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

        fig.add_trace(trace, row=row, col=1)
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
        fig, row = ax

        # Convert linestyle
        dash_map = {
            "-": "solid",
            "--": "dash",
            ":": "dot",
            "-.": "dashdot",
        }
        dash = dash_map.get(linestyle, "solid")

        trace = go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color=color, width=linewidth, dash=dash),
            opacity=alpha,
            name=label or "",
            showlegend=label is not None,
        )

        fig.add_trace(trace, row=row, col=1)
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
        fig, row = ax

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

        fig.add_trace(trace, row=row, col=1)
        return trace

    def axhline(
        self,
        ax: Tuple[go.Figure, int],
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the panel."""
        fig, row = ax

        dash_map = {"-": "solid", "--": "dash", ":": "dot", "-.": "dashdot"}
        dash = dash_map.get(linestyle, "dash")

        fig.add_hline(
            y=y,
            line_dash=dash,
            line_color=color,
            line_width=linewidth,
            row=row,
            col=1,
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
        fig, row = ax

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
            col=1,
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
        fig, row = ax

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
            col=1,
        )

    def set_xlim(self, ax: Tuple[go.Figure, int], left: float, right: float) -> None:
        """Set x-axis limits."""
        fig, row = ax
        xaxis = f"xaxis{row}" if row > 1 else "xaxis"
        fig.update_layout(**{xaxis: dict(range=[left, right])})

    def set_ylim(self, ax: Tuple[go.Figure, int], bottom: float, top: float) -> None:
        """Set y-axis limits."""
        fig, row = ax
        yaxis = f"yaxis{row}" if row > 1 else "yaxis"
        fig.update_layout(**{yaxis: dict(range=[bottom, top])})

    def set_xlabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set x-axis label."""
        fig, row = ax
        xaxis = f"xaxis{row}" if row > 1 else "xaxis"
        fig.update_layout(**{xaxis: dict(title=dict(text=label, font=dict(size=fontsize)))})

    def set_ylabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set y-axis label."""
        fig, row = ax
        yaxis = f"yaxis{row}" if row > 1 else "yaxis"
        fig.update_layout(**{yaxis: dict(title=dict(text=label, font=dict(size=fontsize)))})

    def set_title(
        self, ax: Tuple[go.Figure, int], title: str, fontsize: int = 14
    ) -> None:
        """Set figure title (only works for first panel)."""
        fig, row = ax
        if row == 1:
            fig.update_layout(title=dict(text=title, font=dict(size=fontsize)))

    def create_twin_axis(self, ax: Tuple[go.Figure, int]) -> Tuple[go.Figure, int, str]:
        """Create a secondary y-axis.

        Returns tuple of (figure, row, secondary_yaxis_name).
        """
        fig, row = ax
        secondary_y = f"y{row}2" if row > 1 else "y2"

        # Configure secondary y-axis
        yaxis_name = f"yaxis{row}2" if row > 1 else "yaxis2"
        fig.update_layout(
            **{
                yaxis_name: dict(
                    overlaying=f"y{row}" if row > 1 else "y",
                    side="right",
                )
            }
        )

        return (fig, row, secondary_y)

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

        # Map matplotlib locations to plotly
        loc_map = {
            "upper left": dict(x=0.01, y=0.99, xanchor="left", yanchor="top"),
            "upper right": dict(x=0.99, y=0.99, xanchor="right", yanchor="top"),
            "lower left": dict(x=0.01, y=0.01, xanchor="left", yanchor="bottom"),
            "lower right": dict(x=0.99, y=0.01, xanchor="right", yanchor="bottom"),
        }

        legend_pos = loc_map.get(loc, loc_map["upper left"])
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
        fig, row = ax

        xaxis = f"xaxis{row}" if row > 1 else "xaxis"
        yaxis = f"yaxis{row}" if row > 1 else "yaxis"

        if "top" in spines or "right" in spines:
            # Plotly's template "plotly_white" already hides these
            pass

    def format_xaxis_mb(self, ax: Tuple[go.Figure, int]) -> None:
        """Format x-axis to show megabase values."""
        fig, row = ax
        xaxis = f"xaxis{row}" if row > 1 else "xaxis"

        fig.update_layout(
            **{
                xaxis: dict(
                    tickformat=".2f",
                    ticksuffix=" Mb",
                    tickvals=None,  # Auto
                )
            }
        )

        # Apply custom tick formatting via ticktext/tickvals if needed
        # For now, let plotly auto-format

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

    def finalize_layout(
        self,
        fig: go.Figure,
        left: float = 0.08,
        right: float = 0.95,
        top: float = 0.95,
        bottom: float = 0.1,
        hspace: float = 0.08,
    ) -> None:
        """Adjust layout margins.

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

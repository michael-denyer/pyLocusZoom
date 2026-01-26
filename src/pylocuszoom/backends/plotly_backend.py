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

        # Return (fig, row) tuples for each panel
        # This matches the expected ax parameter format for all methods
        panel_refs = [(fig, row) for row in range(1, n_panels + 1)]
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
            size = max(6, sizes**0.5)  # Approximate conversion
        else:
            size = [max(6, s**0.5) for s in sizes]

        # Build hover template
        if hover_data is not None:
            customdata = hover_data.values
            hover_cols = hover_data.columns.tolist()
            hovertemplate = "<b>%{customdata[0]}</b><br>"
            for i, col in enumerate(hover_cols[1:], 1):
                col_lower = col.lower()
                if (
                    col_lower == "p-value"
                    or col_lower == "pval"
                    or col_lower == "p_value"
                ):
                    hovertemplate += f"{col}: %{{customdata[{i}]:.2e}}<br>"
                elif "r2" in col_lower or "r²" in col_lower or "ld" in col_lower:
                    hovertemplate += f"{col}: %{{customdata[{i}]:.3f}}<br>"
                elif "pos" in col_lower:
                    hovertemplate += f"{col}: %{{customdata[{i}]:,.0f}}<br>"
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
        alpha: float = 1.0,
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
            opacity=alpha,
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
        fig, row = ax

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
        # Convert LaTeX-style labels to Unicode for Plotly
        label = self._convert_label(label)
        fig.update_layout(
            **{xaxis: dict(title=dict(text=label, font=dict(size=fontsize)))}
        )

    def set_ylabel(
        self, ax: Tuple[go.Figure, int], label: str, fontsize: int = 12
    ) -> None:
        """Set y-axis label."""
        fig, row = ax
        yaxis = f"yaxis{row}" if row > 1 else "yaxis"
        # Convert LaTeX-style labels to Unicode for Plotly
        label = self._convert_label(label)
        fig.update_layout(
            **{yaxis: dict(title=dict(text=label, font=dict(size=fontsize)))}
        )

    def _convert_label(self, label: str) -> str:
        """Convert LaTeX-style labels to Unicode for Plotly display."""
        # Common conversions for genomics plots
        conversions = [
            (r"$-\log_{10}$ P", "-log₁₀(P)"),
            (r"$-\log_{10}$", "-log₁₀"),
            (r"\log_{10}", "log₁₀"),
            (r"$r^2$", "r²"),
            (r"$R^2$", "R²"),
        ]
        for latex, unicode_str in conversions:
            if latex in label:
                label = label.replace(latex, unicode_str)
        # Remove any remaining $ markers
        label = label.replace("$", "")
        return label

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
        fig, row = ax

        dash_map = {"-": "solid", "--": "dash", ":": "dot", "-.": "dashdot"}
        dash = dash_map.get(linestyle, "solid")

        trace = go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color=color, width=linewidth, dash=dash),
            opacity=alpha,
            name=label or "",
            showlegend=label is not None,
            yaxis=yaxis_name,
            hoverinfo="skip",
        )

        fig.add_trace(trace, row=row, col=1)
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
        fig, row = ax

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
            yaxis=yaxis_name,
        )

        fig.add_trace(trace, row=row, col=1)
        return trace

    def set_secondary_ylim(
        self,
        ax: Tuple[go.Figure, int],
        bottom: float,
        top: float,
        yaxis_name: str = "y2",
    ) -> None:
        """Set secondary y-axis limits."""
        fig, row = ax
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
        fig, row = ax
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

    def add_ld_legend(
        self,
        ax: Tuple[go.Figure, int],
        ld_bins: List[Tuple[float, str, str]],
        lead_snp_color: str,
    ) -> None:
        """Add LD color legend using invisible scatter traces."""
        fig, row = ax

        # Add LD bin markers (no lead SNP - it's shown in the actual plot)
        for _, label, color in ld_bins:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        symbol="square",
                        size=10,
                        color=color,
                        line=dict(color="black", width=0.5),
                    ),
                    name=label,
                    showlegend=True,
                ),
                row=row,
                col=1,
            )

        # Position legend
        fig.update_layout(
            legend=dict(
                x=0.99,
                y=0.99,
                xanchor="right",
                yanchor="top",
                title=dict(text="r²"),
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        )

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
        # Plotly's template "plotly_white" already hides top/right lines
        # No action needed - method exists for API compatibility
        pass

    def format_xaxis_mb(self, ax: Tuple[go.Figure, int]) -> None:
        """Format x-axis to show megabase values.

        Stores the row for later tick formatting in finalize_layout.
        """
        fig, row = ax
        # Store that this axis needs Mb formatting
        if not hasattr(fig, "_mb_format_rows"):
            fig._mb_format_rows = []
        fig._mb_format_rows.append(row)

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
        """Add eQTL effect size legend using invisible scatter traces."""
        fig, row = ax

        # Positive effects (upward triangles)
        for _, _, label, color in eqtl_positive_bins:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        symbol="triangle-up",
                        size=10,
                        color=color,
                        line=dict(color="black", width=0.5),
                    ),
                    name=label,
                    showlegend=True,
                    legendgroup="eqtl",
                ),
                row=row,
                col=1,
            )

        # Negative effects (downward triangles)
        for _, _, label, color in eqtl_negative_bins:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        symbol="triangle-down",
                        size=10,
                        color=color,
                        line=dict(color="black", width=0.5),
                    ),
                    name=label,
                    showlegend=True,
                    legendgroup="eqtl",
                ),
                row=row,
                col=1,
            )

        # Position legend
        fig.update_layout(
            legend=dict(
                x=0.99,
                y=0.99,
                xanchor="right",
                yanchor="top",
                title=dict(text="eQTL effect"),
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        )

    def add_finemapping_legend(
        self,
        ax: Tuple[go.Figure, int],
        credible_sets: List[int],
        get_color_func: Any,
    ) -> None:
        """Add fine-mapping credible set legend using invisible scatter traces."""
        if not credible_sets:
            return

        fig, row = ax

        for cs_id in credible_sets:
            color = get_color_func(cs_id)
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        symbol="circle",
                        size=10,
                        color=color,
                        line=dict(color="black", width=0.5),
                    ),
                    name=f"CS{cs_id}",
                    showlegend=True,
                    legendgroup="finemapping",
                ),
                row=row,
                col=1,
            )

        # Position legend
        fig.update_layout(
            legend=dict(
                x=0.99,
                y=0.99,
                xanchor="right",
                yanchor="top",
                title=dict(text="Credible sets"),
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        )

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

        loc_map = {
            "upper left": dict(x=0.01, y=0.99, xanchor="left", yanchor="top"),
            "upper right": dict(x=0.99, y=0.99, xanchor="right", yanchor="top"),
            "lower left": dict(x=0.01, y=0.01, xanchor="left", yanchor="bottom"),
            "lower right": dict(x=0.99, y=0.01, xanchor="right", yanchor="bottom"),
        }

        legend_pos = loc_map.get(loc, loc_map["upper right"])
        fig.update_layout(
            legend=dict(
                **legend_pos,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        )

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

            for row in fig._mb_format_rows:
                xaxis_name = f"xaxis{row}" if row > 1 else "xaxis"
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
                    # Create nice tick values in Mb
                    x_min_mb = x_range[0] / 1e6
                    x_max_mb = x_range[1] / 1e6
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

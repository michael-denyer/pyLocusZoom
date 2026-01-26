"""Bokeh backend for pyLocusZoom.

Interactive backend with hover tooltips, well-suited for dashboards.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
from bokeh.io import export_png, export_svgs, output_file, save, show
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, HoverTool, Legend, LegendItem, Span
from bokeh.plotting import figure


class BokehBackend:
    """Bokeh backend for interactive plot generation.

    Produces interactive HTML plots suitable for embedding in web
    applications and dashboards.
    """

    def __init__(self) -> None:
        """Initialize the bokeh backend."""
        self._marker_map = {
            "o": "circle",
            "D": "diamond",
            "s": "square",
            "^": "triangle",
            "v": "inverted_triangle",
        }

    def create_figure(
        self,
        n_panels: int,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[Any, List[figure]]:
        """Create a layout with multiple panels.

        Args:
            n_panels: Number of vertical panels.
            height_ratios: Relative heights for each panel.
            figsize: Figure size as (width, height) in inches.
            sharex: Whether panels share the x-axis.

        Returns:
            Tuple of (layout, list of figure objects).
        """
        # Convert inches to pixels
        width_px = int(figsize[0] * 100)
        total_height = int(figsize[1] * 100)

        # Calculate individual heights
        total_ratio = sum(height_ratios)
        heights = [int(total_height * r / total_ratio) for r in height_ratios]

        figures = []
        x_range = None

        for i, h in enumerate(heights):
            p = figure(
                width=width_px,
                height=h,
                x_range=x_range if sharex and x_range else None,
                tools="pan,wheel_zoom,box_zoom,reset,save",
                toolbar_location="above" if i == 0 else None,
            )

            # Share x_range for subsequent figures
            if sharex and x_range is None:
                x_range = p.x_range

            # Style
            p.grid.grid_line_alpha = 0.3
            p.outline_line_color = None

            figures.append(p)

        # Create column layout
        layout = column(*figures, sizing_mode="fixed")

        return layout, figures

    def scatter(
        self,
        ax: figure,
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
        """Create a scatter plot on the given figure."""
        # Prepare data source
        data = {"x": x.values, "y": y.values}

        # Handle colors
        if isinstance(colors, str):
            data["color"] = [colors] * len(x)
        else:
            data["color"] = list(colors) if hasattr(colors, "tolist") else colors

        # Handle sizes (convert from area to diameter)
        if isinstance(sizes, (int, float)):
            bokeh_size = max(6, sizes ** 0.5)
            data["size"] = [bokeh_size] * len(x)
        else:
            data["size"] = [max(6, s ** 0.5) for s in sizes]

        # Add hover data
        tooltips = []
        if hover_data is not None:
            for col in hover_data.columns:
                data[col] = hover_data[col].values
                if "p" in col.lower():
                    tooltips.append((col, "@{" + col + "}{0.2e}"))
                elif "r2" in col.lower() or "ld" in col.lower():
                    tooltips.append((col, "@{" + col + "}{0.3f}"))
                else:
                    tooltips.append((col, f"@{col}"))

        source = ColumnDataSource(data)

        # Get marker type
        marker_type = self._marker_map.get(marker, "circle")

        # Create scatter
        renderer = getattr(ax, marker_type)(
            "x",
            "y",
            source=source,
            size="size",
            fill_color="color",
            line_color=edgecolor,
            line_width=linewidth,
            legend_label=label if label else None,
        )

        # Add hover tool if we have hover data
        if tooltips:
            hover = HoverTool(
                tooltips=tooltips,
                renderers=[renderer],
                mode="mouse",
            )
            ax.add_tools(hover)

        return renderer

    def line(
        self,
        ax: figure,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        zorder: int = 1,
        label: Optional[str] = None,
    ) -> Any:
        """Create a line plot on the given figure."""
        # Convert linestyle
        dash_map = {
            "-": "solid",
            "--": "dashed",
            ":": "dotted",
            "-.": "dashdot",
        }
        line_dash = dash_map.get(linestyle, "solid")

        return ax.line(
            x.values,
            y.values,
            line_color=color,
            line_width=linewidth,
            line_alpha=alpha,
            line_dash=line_dash,
            legend_label=label if label else None,
        )

    def fill_between(
        self,
        ax: figure,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> Any:
        """Fill area between two y-values."""
        # Convert to arrays
        x_arr = x.values
        if isinstance(y1, (int, float)):
            y1_arr = [y1] * len(x_arr)
        else:
            y1_arr = y1.values if hasattr(y1, "values") else list(y1)

        if isinstance(y2, (int, float)):
            y2_arr = [y2] * len(x_arr)
        else:
            y2_arr = y2.values if hasattr(y2, "values") else list(y2)

        return ax.varea(
            x=x_arr,
            y1=y1_arr,
            y2=y2_arr,
            fill_color=color,
            fill_alpha=alpha,
        )

    def axhline(
        self,
        ax: figure,
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the figure."""
        dash_map = {"-": "solid", "--": "dashed", ":": "dotted", "-.": "dashdot"}
        line_dash = dash_map.get(linestyle, "dashed")

        span = Span(
            location=y,
            dimension="width",
            line_color=color,
            line_dash=line_dash,
            line_width=linewidth,
        )
        ax.add_layout(span)
        return span

    def add_text(
        self,
        ax: figure,
        x: float,
        y: float,
        text: str,
        fontsize: int = 10,
        ha: str = "center",
        va: str = "bottom",
        rotation: float = 0,
        color: str = "black",
    ) -> Any:
        """Add text annotation to figure."""
        from bokeh.models import Label

        # Map alignment
        anchor_map = {
            ("center", "bottom"): ("center", "bottom"),
            ("center", "top"): ("center", "top"),
            ("left", "bottom"): ("left", "bottom"),
            ("right", "bottom"): ("right", "bottom"),
        }
        text_align, text_baseline = anchor_map.get((ha, va), ("center", "bottom"))

        label = Label(
            x=x,
            y=y,
            text=text,
            text_font_size=f"{fontsize}pt",
            text_color=color,
            text_align=text_align,
            text_baseline=text_baseline,
            angle=rotation,
            angle_units="deg",
        )
        ax.add_layout(label)
        return label

    def add_rectangle(
        self,
        ax: figure,
        xy: Tuple[float, float],
        width: float,
        height: float,
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a rectangle to the figure."""
        from bokeh.models import Rect

        x_center = xy[0] + width / 2
        y_center = xy[1] + height / 2

        return ax.rect(
            x=[x_center],
            y=[y_center],
            width=[width],
            height=[height],
            fill_color=facecolor,
            line_color=edgecolor,
            line_width=linewidth,
        )

    def set_xlim(self, ax: figure, left: float, right: float) -> None:
        """Set x-axis limits."""
        ax.x_range.start = left
        ax.x_range.end = right

    def set_ylim(self, ax: figure, bottom: float, top: float) -> None:
        """Set y-axis limits."""
        ax.y_range.start = bottom
        ax.y_range.end = top

    def set_xlabel(self, ax: figure, label: str, fontsize: int = 12) -> None:
        """Set x-axis label."""
        ax.xaxis.axis_label = label
        ax.xaxis.axis_label_text_font_size = f"{fontsize}pt"

    def set_ylabel(self, ax: figure, label: str, fontsize: int = 12) -> None:
        """Set y-axis label."""
        ax.yaxis.axis_label = label
        ax.yaxis.axis_label_text_font_size = f"{fontsize}pt"

    def set_title(self, ax: figure, title: str, fontsize: int = 14) -> None:
        """Set figure title."""
        ax.title.text = title
        ax.title.text_font_size = f"{fontsize}pt"

    def create_twin_axis(self, ax: figure) -> Any:
        """Create a secondary y-axis.

        Returns a dict with configuration for extra_y_ranges.
        """
        from bokeh.models import LinearAxis, Range1d

        # Add a second y-axis
        ax.extra_y_ranges = {"secondary": Range1d(start=0, end=100)}
        ax.add_layout(LinearAxis(y_range_name="secondary"), "right")

        return "secondary"

    def add_legend(
        self,
        ax: figure,
        handles: List[Any],
        labels: List[str],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> Any:
        """Configure legend on the figure."""
        # Bokeh handles legend automatically from legend_label
        # Just configure position

        loc_map = {
            "upper left": "top_left",
            "upper right": "top_right",
            "lower left": "bottom_left",
            "lower right": "bottom_right",
        }

        ax.legend.location = loc_map.get(loc, "top_left")
        if title:
            ax.legend.title = title
        ax.legend.background_fill_alpha = 0.9
        ax.legend.border_line_color = "black"

        return ax.legend

    def hide_spines(self, ax: figure, spines: List[str]) -> None:
        """Hide specified axis spines."""
        # Bokeh doesn't have spines in the same way
        # We can hide axis lines
        if "top" in spines:
            ax.xaxis.visible = ax.xaxis.visible  # Keep visible but could customize
        if "right" in spines:
            ax.yaxis.visible = ax.yaxis.visible

    def format_xaxis_mb(self, ax: figure) -> None:
        """Format x-axis to show megabase values."""
        from bokeh.models import NumeralTickFormatter

        ax.xaxis.formatter = NumeralTickFormatter(format="0.00")
        ax.xaxis.axis_label = ax.xaxis.axis_label or "Position (Mb)"

        # We need to scale values or use a custom formatter
        # For now, assume values are already in bp and need /1e6
        from bokeh.models import FuncTickFormatter

        ax.xaxis.formatter = FuncTickFormatter(
            code="return (tick / 1e6).toFixed(2);"
        )

    def save(
        self,
        fig: Any,
        path: str,
        dpi: int = 150,
        bbox_inches: str = "tight",
    ) -> None:
        """Save figure to file.

        Supports .html for interactive and .png for static.
        """
        if path.endswith(".html"):
            output_file(path)
            save(fig)
        elif path.endswith(".png"):
            export_png(fig, filename=path)
        elif path.endswith(".svg"):
            export_svgs(fig, filename=path)
        else:
            # Default to HTML
            output_file(path)
            save(fig)

    def show(self, fig: Any) -> None:
        """Display the figure."""
        show(fig)

    def close(self, fig: Any) -> None:
        """Close the figure (no-op for bokeh)."""
        pass

    def finalize_layout(
        self,
        fig: Any,
        left: float = 0.08,
        right: float = 0.95,
        top: float = 0.95,
        bottom: float = 0.1,
        hspace: float = 0.08,
    ) -> None:
        """Adjust layout (limited support in bokeh).

        Bokeh layouts are mostly automatic.
        """
        # Bokeh handles layout differently - column spacing is fixed
        pass

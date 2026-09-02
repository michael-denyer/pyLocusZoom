"""Bokeh backend for pyLocusZoom.

Interactive backend with hover tooltips, well-suited for dashboards.
"""

import math
from typing import Any, List, NamedTuple, Optional, Tuple, Union

import numpy as np
import pandas as pd
from bokeh.layouts import column, row
from bokeh.models import (
    BasicTicker,
    BoxAnnotation,
    ColorBar,
    ColumnDataSource,
    CustomJSTickFormatter,
    DataRange1d,
    HoverTool,
    Label,
    Legend,
    LegendItem,
    LinearAxis,
    LinearColorMapper,
    Range1d,
    Span,
    Whisker,
)
from bokeh.models.layouts import Column
from bokeh.plotting import figure
from matplotlib.colors import LinearSegmentedColormap, to_hex

from . import convert_latex_to_unicode, register_backend
from ._coerce import (
    broadcast,
    marker_colors,
    marker_diameter,
    per_point,
    pixels,
    split_pixels,
)
from .composition import LegendEntry
from .hover import bokeh_tooltips

# Style mappings (matplotlib -> Bokeh)
_MARKER_MAP = {
    "o": "circle",
    "D": "diamond",
    "s": "square",
    "^": "triangle",
    "v": "inverted_triangle",
}
# Matplotlib legend `loc` vocabulary mapped to Bokeh legend locations.
# "best" has no Bokeh equivalent, so it takes the upper-right default.
_LEGEND_LOCATIONS = {
    "best": "top_right",
    "upper right": "top_right",
    "upper left": "top_left",
    "upper center": "top_center",
    "lower right": "bottom_right",
    "lower left": "bottom_left",
    "lower center": "bottom_center",
    "center right": "center_right",
    "center left": "center_left",
    "right": "center_right",
    "center": "center",
}
_DASH_MAP = {
    "-": "solid",
    "--": "dashed",
    ":": "dotted",
    "-.": "dashdot",
}
# matplotlib va vocabulary to bokeh text_baseline; "top" and "bottom" pass through.
_BASELINE_MAP = {"center": "middle", "baseline": "alphabetic"}
# Namespaces hover columns in a ColumnDataSource so a hover column named "x"
# or "size" cannot shadow the keys scatter() sets for geometry and styling.
_HOVER_KEY_PREFIX = "hover_"
# Names of the y-ranges a glyph can be drawn against.
_DEFAULT_RANGE = "default"
_SECONDARY_RANGE = "secondary"


class _SecondaryAxis(NamedTuple):
    """A figure's secondary y-axis, as returned by ``create_twin_axis``.

    Carries the axis and the range ``create_twin_axis`` added, so labelling and
    scaling reach them directly instead of searching the figure for them.
    """

    figure: figure
    axis: LinearAxis
    y_range: Range1d
    name: str = _SECONDARY_RANGE


def _draw_target(ax: Union[figure, _SecondaryAxis]) -> Tuple[figure, str]:
    """The figure to draw on and the y-range to draw against."""
    if isinstance(ax, _SecondaryAxis):
        return ax.figure, ax.name
    return ax, _DEFAULT_RANGE


@register_backend("bokeh")
class BokehBackend:
    """Bokeh backend for interactive plot generation.

    Produces interactive HTML plots suitable for embedding in web
    applications and dashboards.
    """

    @property
    def supports_hover(self) -> bool:
        """Bokeh supports hover tooltips."""
        return True

    def create_figure(
        self,
        n_panels: int,
        height_ratios: List[float],
        figsize: Tuple[float, float],
        sharex: bool = True,
    ) -> Tuple[Any, List[figure]]:
        """Create a layout with multiple panels."""
        width_px, total_height = pixels(figsize)
        heights = split_pixels(total_height, height_ratios, len(height_ratios))

        figures = []
        x_range = DataRange1d() if sharex else None

        for i, h in enumerate(heights):
            p = figure(
                width=width_px,
                height=h,
                x_range=x_range if sharex else DataRange1d(),
                tools="pan,wheel_zoom,box_zoom,reset,save",
                toolbar_location="above" if i == 0 else None,
            )

            # Style - no grid lines, black axes for clean LocusZoom appearance
            p.grid.visible = False
            p.outline_line_color = None
            p.xaxis.axis_line_color = "black"
            p.yaxis.axis_line_color = "black"
            p.xaxis.minor_tick_line_color = None
            p.yaxis.minor_tick_line_color = None

            figures.append(p)

        # Create column layout (use default sizing mode to avoid validation warnings)
        layout = column(*figures)

        return layout, figures

    def create_figure_grid(
        self,
        n_rows: int,
        n_cols: int,
        width_ratios: Optional[List[float]] = None,
        height_ratios: Optional[List[float]] = None,
        figsize: Tuple[float, float] = (12.0, 8.0),
    ) -> Tuple[Any, List[figure]]:
        """Create a layout with a grid of subplots."""
        width_px, height_px = pixels(figsize)
        widths = split_pixels(width_px, width_ratios, n_cols)
        heights = split_pixels(height_px, height_ratios, n_rows)

        figures = []
        rows = []

        for i in range(n_rows):
            row_figures = []
            for j in range(n_cols):
                p = figure(
                    width=widths[j],
                    height=heights[i],
                    tools="pan,wheel_zoom,box_zoom,reset,save",
                    toolbar_location="above" if i == 0 and j == 0 else None,
                )

                # Style
                p.grid.visible = False
                p.outline_line_color = None
                p.xaxis.axis_line_color = "black"
                p.yaxis.axis_line_color = "black"
                p.xaxis.minor_tick_line_color = None
                p.yaxis.minor_tick_line_color = None

                row_figures.append(p)
                figures.append(p)

            rows.append(row(*row_figures))

        layout = column(*rows)

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
    ) -> Any:
        """Create a scatter plot on the given figure."""
        # Prepare data source
        data = {"x": x.values, "y": y.values}

        data["color"] = per_point(marker_colors(colors), len(x))
        data["size"] = per_point(marker_diameter(sizes), len(x))

        # Add hover data with namespaced keys to avoid collisions
        # with internal keys (x, y, color, size)
        tooltips = []
        if hover_data is not None:
            for col in hover_data.columns:
                data[f"{_HOVER_KEY_PREFIX}{col}"] = hover_data[col].values
            tooltips = bokeh_tooltips(hover_data, key_prefix=_HOVER_KEY_PREFIX)

        source = ColumnDataSource(data)

        marker_type = _MARKER_MAP.get(marker, "circle")

        # Create scatter using scatter() method (Bokeh 3.4+ preferred API)
        renderer = ax.scatter(
            "x",
            "y",
            source=source,
            marker=marker_type,
            size="size",
            fill_color="color",
            line_color=edgecolor,
            line_width=linewidth,
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
        ax: Union[figure, _SecondaryAxis],
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        zorder: int = 1,
    ) -> Any:
        """Create a line plot on the given figure or secondary axis."""
        target, y_range_name = _draw_target(ax)

        return target.line(
            x.values,
            y.values,
            line_color=color,
            line_width=linewidth,
            line_alpha=alpha,
            line_dash=_DASH_MAP.get(linestyle, "solid"),
            y_range_name=y_range_name,
        )

    def fill_between(
        self,
        ax: Union[figure, _SecondaryAxis],
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        zorder: int = 0,
    ) -> Any:
        """Fill area between two y-values."""
        target, y_range_name = _draw_target(ax)
        x_arr = x.values
        return target.varea(
            x=x_arr,
            y1=broadcast(y1, len(x_arr)),
            y2=broadcast(y2, len(x_arr)),
            fill_color=color,
            fill_alpha=alpha,
            y_range_name=y_range_name,
        )

    def axhline(
        self,
        ax: figure,
        y: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a horizontal line across the figure."""
        line_dash = _DASH_MAP.get(linestyle, "dashed")

        span = Span(
            location=y,
            dimension="width",
            line_color=color,
            line_dash=line_dash,
            line_width=linewidth,
            line_alpha=alpha,
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
        label = Label(
            x=x,
            y=y,
            text=text,
            text_font_size=f"{fontsize}pt",
            text_color=color,
            text_align=ha,
            text_baseline=_BASELINE_MAP.get(va, va),
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
        facecolor: Optional[str] = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a rectangle to the figure."""
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

    def add_polygon(
        self,
        ax: figure,
        points: List[List[float]],
        facecolor: str = "blue",
        edgecolor: str = "black",
        linewidth: float = 0.5,
        zorder: int = 2,
    ) -> Any:
        """Add a polygon (e.g., triangle for strand arrows) to the figure."""
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]

        # Bokeh patch() uses x/y (singular) for single polygon
        return ax.patch(
            x=xs,
            y=ys,
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
        label = convert_latex_to_unicode(label)
        ax.xaxis.axis_label = label
        ax.xaxis.axis_label_text_font_size = f"{fontsize}pt"

    def set_ylabel(self, ax: figure, label: str, fontsize: int = 12) -> None:
        """Set y-axis label."""
        label = convert_latex_to_unicode(label)
        ax.yaxis.axis_label = label
        ax.yaxis.axis_label_text_font_size = f"{fontsize}pt"

    def set_yticks(
        self,
        ax: figure,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
    ) -> None:
        """Set y-axis tick positions and labels."""
        ax.yaxis.ticker = positions
        ax.yaxis.major_label_overrides = {
            pos: label for pos, label in zip(positions, labels)
        }
        ax.yaxis.major_label_text_font_size = f"{fontsize}pt"

    def set_xticks(
        self,
        ax: figure,
        positions: List[float],
        labels: List[str],
        fontsize: int = 10,
        rotation: int = 0,
        ha: str = "center",
    ) -> None:
        """Set x-axis tick positions and labels."""
        ax.xaxis.ticker = positions
        ax.xaxis.major_label_overrides = {
            pos: label for pos, label in zip(positions, labels)
        }
        ax.xaxis.major_label_text_font_size = f"{fontsize}pt"
        if rotation:
            ax.xaxis.major_label_orientation = math.radians(rotation)

    def set_title(self, ax: figure, title: str, fontsize: int = 14) -> None:
        """Set figure title."""
        ax.title.text = title
        ax.title.text_font_size = f"{fontsize}pt"

    def set_suptitle(self, fig: Any, title: str, fontsize: int = 14) -> None:
        """Set overall figure title.

        For Bokeh layouts, add title to the first figure in the layout.
        """
        if isinstance(fig, Column) and len(fig.children) > 0:
            first_child = fig.children[0]
            if hasattr(first_child, "title"):
                first_child.title.text = title
                first_child.title.text_font_size = f"{fontsize}pt"
        elif hasattr(fig, "title"):
            fig.title.text = title
            fig.title.text_font_size = f"{fontsize}pt"

    def create_twin_axis(self, ax: figure) -> _SecondaryAxis:
        """Create a secondary y-axis and return its handle."""
        # Add a second y-axis without tick marks (cleaner look)
        y_range = Range1d(start=0, end=100)
        ax.extra_y_ranges = {_SECONDARY_RANGE: y_range}
        secondary_axis = LinearAxis(
            y_range_name=_SECONDARY_RANGE,
            major_tick_line_color=None,  # Hide major ticks
            minor_tick_line_color=None,  # Hide minor ticks
            major_label_text_font_size="0pt",  # Hide tick labels
        )
        ax.add_layout(secondary_axis, "right")

        return _SecondaryAxis(ax, secondary_axis, y_range)

    def set_secondary_ylim(
        self,
        secondary: _SecondaryAxis,
        bottom: float,
        top: float,
    ) -> None:
        """Set secondary y-axis limits."""
        secondary.y_range.start = bottom
        secondary.y_range.end = top

    def set_secondary_ylabel(
        self,
        secondary: _SecondaryAxis,
        label: str,
        color: str = "black",
        fontsize: int = 10,
    ) -> None:
        """Set secondary y-axis label."""
        secondary.axis.axis_label = convert_latex_to_unicode(label)
        secondary.axis.axis_label_text_font_size = f"{fontsize}pt"
        secondary.axis.axis_label_text_color = color
        secondary.axis.major_label_text_color = color

    def add_panel_label(
        self,
        ax: figure,
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        # Use screen coordinates so the label works regardless of whether
        # the data range has been resolved (DataRange1d starts as None)
        x_px = int(x_frac * ax.width)
        y_px = int(y_frac * ax.height)

        label_obj = Label(
            x=x_px,
            y=y_px,
            x_units="screen",
            y_units="screen",
            text=label,
            text_font_size="12px",
            text_font_style="bold",
        )
        ax.add_layout(label_obj)

    def _ensure_legend_range(self, ax: figure) -> Any:
        """Ensure legend range exists and return a dummy data source.

        Creates a separate y-range for legend glyphs so they don't affect
        the main plot's axis scaling.
        """
        if "legend_range" not in ax.extra_y_ranges:
            ax.extra_y_ranges["legend_range"] = Range1d(start=0, end=1)
        return ColumnDataSource(data={"x": [0], "y": [0]})

    def _add_legend_item(
        self,
        ax: figure,
        source: Any,
        label: str,
        color: str,
        marker: str,
        size: int = 14,
        edgecolor: str = "black",
    ) -> Any:
        """Create an invisible scatter renderer for a legend entry."""
        renderer = ax.scatter(
            x="x",
            y="y",
            source=source,
            marker=marker,
            size=size,
            fill_color=color,
            line_color=edgecolor,
            line_width=0.5,
            y_range_name="legend_range",
            visible=False,
        )
        return LegendItem(label=label, renderers=[renderer])

    def _create_legend(
        self, ax: figure, items: List[Any], title: str, loc: str
    ) -> None:
        """Create and add a styled legend to the figure."""
        legend = Legend(
            items=items,
            location=_LEGEND_LOCATIONS.get(loc, "top_right"),
            title=convert_latex_to_unicode(title),
            background_fill_alpha=0.9,
            border_line_color="black",
            spacing=0,
            padding=4,
            label_height=12,
            glyph_height=12,
        )
        ax.add_layout(legend)

    def add_legend(
        self,
        ax: figure,
        entries: List[LegendEntry],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> None:
        """Render legend entries as a Bokeh legend using invisible glyphs."""
        source = self._ensure_legend_range(ax)
        items = []
        for entry in entries:
            marker = (
                "square"
                if entry.marker == "patch"
                else _MARKER_MAP.get(entry.marker, "circle")
            )
            items.append(
                self._add_legend_item(
                    ax,
                    source,
                    entry.label,
                    entry.color,
                    marker,
                    edgecolor=entry.edgecolor or "black",
                )
            )
        self._create_legend(ax, items, title or "", loc)

    def hide_yaxis(self, ax: figure) -> None:
        """Hide y-axis ticks, labels, line, and grid for gene track panels."""
        ax.yaxis.visible = False
        ax.ygrid.visible = False

    def format_xaxis_mb(self, ax: figure) -> None:
        """Format x-axis to show megabase values."""
        ax.xaxis.formatter = CustomJSTickFormatter(
            code="return (tick / 1e6).toFixed(2);"
        )

    def axvline(
        self,
        ax: figure,
        x: float,
        color: str = "grey",
        linestyle: str = "--",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        zorder: int = 1,
    ) -> Any:
        """Add a vertical line across the figure."""
        line_dash = _DASH_MAP.get(linestyle, "dashed")

        span = Span(
            location=x,
            dimension="height",
            line_color=color,
            line_dash=line_dash,
            line_width=linewidth,
            line_alpha=alpha,
        )
        ax.add_layout(span)
        return span

    def hbar(
        self,
        ax: figure,
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
        # Convert left to array if scalar
        if isinstance(left, (int, float)):
            left_arr = [left] * len(y)
        else:
            left_arr = list(left) if hasattr(left, "tolist") else left

        # Calculate right edge
        right_arr = [left_val + w for left_val, w in zip(left_arr, width)]

        return ax.hbar(
            y=y.values,
            right=right_arr,
            left=left_arr,
            height=height,
            fill_color=color,
            line_color=edgecolor,
            line_width=linewidth,
        )

    def errorbar_h(
        self,
        ax: figure,
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
        # Calculate bounds
        lower = x - xerr_lower
        upper = x + xerr_upper

        source = ColumnDataSource(
            data={
                "y": y.values,
                "lower": lower.values,
                "upper": upper.values,
            }
        )

        # Add horizontal whisker
        whisker = Whisker(
            source=source,
            base="y",
            lower="lower",
            upper="upper",
            dimension="width",
            line_color=color,
            line_width=linewidth,
        )
        ax.add_layout(whisker)
        return whisker

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

    def add_region_highlight(
        self,
        fig: Any,
        axes: List[figure],
        x_start: float,
        x_end: float,
        color: str = "yellow",
        alpha: float = 0.3,
    ) -> None:
        """Highlight an x-range across multiple Bokeh panels."""
        for ax in axes:
            ax.add_layout(
                BoxAnnotation(
                    left=x_start,
                    right=x_end,
                    fill_color=color,
                    fill_alpha=alpha,
                )
            )

    def add_heatmap(
        self,
        ax: figure,
        data: Any,
        x_coords: List[float],
        y_coords: List[float],
        cmap_colors: Optional[List[str]] = None,
        vmin: float = 0.0,
        vmax: float = 1.0,
    ) -> Any:
        """Render a heatmap of an already-shaped matrix."""
        if cmap_colors is None:
            cmap_colors = ["#FFFFFF", "#FF0000"]

        # Create custom palette from start to end color
        # For a simple 2-color gradient, create a palette of intermediate colors
        palette = _create_color_palette(cmap_colors[0], cmap_colors[1], 256)

        mapper = LinearColorMapper(
            palette=palette,
            low=vmin,
            high=vmax,
            nan_color="#808080",  # Grey for missing
        )

        # A masked cell is one the caller shaped out, so no rect is emitted for
        # it. An unmasked NaN is missing data and still draws, in nan_color.
        masked = np.ma.getmaskarray(np.ma.asarray(data))
        n = data.shape[0]
        cells = [(i, j) for i in range(n) for j in range(n) if not masked[i, j]]
        xs = [x_coords[j] for _, j in cells]
        ys = [y_coords[i] for i, _ in cells]
        values = [float(data[i, j]) for i, j in cells]

        # Compute per-cell widths and heights based on actual coordinate spacing.
        # Uses midpoints between adjacent coordinates for cell boundaries.
        def _cell_sizes(coords: List[float]) -> List[float]:
            if len(coords) <= 1:
                return [1.0] * len(coords)
            sizes = []
            for i in range(len(coords)):
                left = (
                    (coords[i - 1] + coords[i]) / 2
                    if i > 0
                    else coords[i] - (coords[1] - coords[0]) / 2
                )
                right = (
                    (coords[i] + coords[i + 1]) / 2
                    if i < len(coords) - 1
                    else coords[i] + (coords[-1] - coords[-2]) / 2
                )
                sizes.append(abs(right - left))
            return sizes

        x_sizes = _cell_sizes(x_coords)
        y_sizes = _cell_sizes(y_coords)

        widths = [x_sizes[j] for _, j in cells]
        heights = [y_sizes[i] for i, _ in cells]

        source = ColumnDataSource(
            {"x": xs, "y": ys, "value": values, "w": widths, "h": heights}
        )

        ax.rect(
            x="x",
            y="y",
            width="w",
            height="h",
            fill_color={"field": "value", "transform": mapper},
            line_color=None,
            source=source,
        )
        return mapper

    def add_colorbar(
        self,
        ax: figure,
        mappable: Any,
        label: str = "R²",
        orientation: str = "vertical",
    ) -> Any:
        """Add colorbar legend for heatmap."""
        color_bar = ColorBar(
            color_mapper=mappable,
            ticker=BasicTicker(),
            label_standoff=6,
            title=label,
            orientation=orientation,
        )
        ax.add_layout(color_bar, "right")
        return color_bar


def _create_color_palette(start_color: str, end_color: str, n_colors: int) -> List[str]:
    """Create a linear color palette between two colors.

    Args:
        start_color: Starting color (hex or named CSS color).
        end_color: Ending color (hex or named CSS color).
        n_colors: Number of colors in palette.

    Returns:
        List of hex color strings.
    """
    cmap = LinearSegmentedColormap.from_list("ld", [start_color, end_color], N=n_colors)
    return [to_hex(cmap(i)) for i in range(n_colors)]

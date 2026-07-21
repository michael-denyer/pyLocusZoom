"""Bokeh backend for pyLocusZoom.

Interactive backend with hover tooltips, well-suited for dashboards.
"""

import math
from typing import Any, List, Optional, Tuple, Union

import pandas as pd
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, DataRange1d, HoverTool, Span
from bokeh.plotting import figure

from . import convert_latex_to_unicode, register_backend
from .composition import LegendEntry

# Style mappings (matplotlib -> Bokeh)
_MARKER_MAP = {
    "o": "circle",
    "D": "diamond",
    "s": "square",
    "^": "triangle",
    "v": "inverted_triangle",
}
_DASH_MAP = {
    "-": "solid",
    "--": "dashed",
    ":": "dotted",
    "-.": "dashdot",
}


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
        """Create a layout with a grid of subplots.

        Args:
            n_rows: Number of rows.
            n_cols: Number of columns.
            width_ratios: Relative widths for columns.
            height_ratios: Relative heights for rows.
            figsize: Figure size as (width, height).

        Returns:
            Tuple of (layout, flattened list of figure objects).
        """
        width_px = int(figsize[0] * 100)
        height_px = int(figsize[1] * 100)

        # Calculate widths
        if width_ratios is not None:
            total_w = sum(width_ratios)
            widths = [int(width_px * w / total_w) for w in width_ratios]
        else:
            widths = [width_px // n_cols] * n_cols

        # Calculate heights
        if height_ratios is not None:
            total_h = sum(height_ratios)
            heights = [int(height_px * h / total_h) for h in height_ratios]
        else:
            heights = [height_px // n_rows] * n_rows

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
            bokeh_size = max(6, sizes**0.5)
            data["size"] = [bokeh_size] * len(x)
        else:
            data["size"] = [max(6, s**0.5) for s in sizes]

        # Add hover data with namespaced keys to avoid collisions
        # with internal keys (x, y, color, size)
        tooltips = []
        if hover_data is not None:
            for col in hover_data.columns:
                safe_key = f"hover_{col}"
                data[safe_key] = hover_data[col].values
                col_lower = col.lower()
                if col_lower in ("p-value", "pval", "p_value"):
                    tooltips.append((col, "@{" + safe_key + "}{0.2e}"))
                elif any(kw in col_lower for kw in ("r2", "r²", "ld")):
                    tooltips.append((col, "@{" + safe_key + "}{0.3f}"))
                elif "pos" in col_lower:
                    tooltips.append((col, "@{" + safe_key + "}{0,0}"))
                else:
                    tooltips.append((col, f"@{{{safe_key}}}"))

        source = ColumnDataSource(data)

        marker_type = _MARKER_MAP.get(marker, "circle")

        # Create scatter using scatter() method (Bokeh 3.4+ preferred API)
        scatter_kwargs = {
            "source": source,
            "marker": marker_type,
            "size": "size",
            "fill_color": "color",
            "line_color": edgecolor,
            "line_width": linewidth,
        }
        if label:
            scatter_kwargs["legend_label"] = label

        renderer = ax.scatter("x", "y", **scatter_kwargs)

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
        line_dash = _DASH_MAP.get(linestyle, "solid")

        line_kwargs = {
            "line_color": color,
            "line_width": linewidth,
            "line_alpha": alpha,
            "line_dash": line_dash,
        }
        if label:
            line_kwargs["legend_label"] = label

        return ax.line(x.values, y.values, **line_kwargs)

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
        label = self._convert_label(label)
        ax.xaxis.axis_label = label
        ax.xaxis.axis_label_text_font_size = f"{fontsize}pt"

    def set_ylabel(self, ax: figure, label: str, fontsize: int = 12) -> None:
        """Set y-axis label."""
        label = self._convert_label(label)
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

    def _convert_label(self, label: str) -> str:
        """Convert LaTeX-style labels to Unicode for Bokeh display."""
        return convert_latex_to_unicode(label)

    def set_title(self, ax: figure, title: str, fontsize: int = 14) -> None:
        """Set figure title."""
        ax.title.text = title
        ax.title.text_font_size = f"{fontsize}pt"

    def set_suptitle(self, fig: Any, title: str, fontsize: int = 14) -> None:
        """Set overall figure title.

        For Bokeh layouts, add title to the first figure in the layout.
        """
        from bokeh.models.layouts import Column

        if isinstance(fig, Column) and len(fig.children) > 0:
            first_child = fig.children[0]
            if hasattr(first_child, "title"):
                first_child.title.text = title
                first_child.title.text_font_size = f"{fontsize}pt"
        elif hasattr(fig, "title"):
            fig.title.text = title
            fig.title.text_font_size = f"{fontsize}pt"

    def create_twin_axis(self, ax: figure) -> Any:
        """Create a secondary y-axis.

        Returns an opaque ``(ax, yaxis_name)`` handle for the ``*_secondary``
        primitives.
        """
        from bokeh.models import LinearAxis, Range1d

        # Add a second y-axis without tick marks (cleaner look)
        ax.extra_y_ranges = {"secondary": Range1d(start=0, end=100)}
        secondary_axis = LinearAxis(
            y_range_name="secondary",
            major_tick_line_color=None,  # Hide major ticks
            minor_tick_line_color=None,  # Hide minor ticks
            major_label_text_font_size="0pt",  # Hide tick labels
        )
        ax.add_layout(secondary_axis, "right")

        return (ax, "secondary")

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
        line_dash = _DASH_MAP.get(linestyle, "solid")

        return ax.line(
            x.values,
            y.values,
            line_color=color,
            line_width=linewidth,
            line_alpha=alpha,
            line_dash=line_dash,
            y_range_name=yaxis_name,
        )

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
            y_range_name=yaxis_name,
        )

    def set_secondary_ylim(
        self,
        secondary: Any,
        bottom: float,
        top: float,
    ) -> None:
        """Set secondary y-axis limits."""
        ax, yaxis_name = secondary
        if yaxis_name in ax.extra_y_ranges:
            ax.extra_y_ranges[yaxis_name].start = bottom
            ax.extra_y_ranges[yaxis_name].end = top

    def set_secondary_ylabel(
        self,
        secondary: Any,
        label: str,
        color: str = "black",
        fontsize: int = 10,
    ) -> None:
        """Set secondary y-axis label."""
        ax, yaxis_name = secondary
        label = self._convert_label(label)
        # Find the secondary axis and update its label
        for renderer in ax.right:
            if (
                hasattr(renderer, "y_range_name")
                and renderer.y_range_name == yaxis_name
            ):
                renderer.axis_label = label
                renderer.axis_label_text_font_size = f"{fontsize}pt"
                renderer.axis_label_text_color = color
                renderer.major_label_text_color = color
                break

    def add_panel_label(
        self,
        ax: figure,
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        from bokeh.models import Label

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
        from bokeh.models import Range1d

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
    ) -> Any:
        """Create an invisible scatter renderer for a legend entry."""
        from bokeh.models import LegendItem

        renderer = ax.scatter(
            x="x",
            y="y",
            source=source,
            marker=marker,
            size=size,
            fill_color=color,
            line_color="black",
            line_width=0.5,
            y_range_name="legend_range",
            visible=False,
        )
        return LegendItem(label=label, renderers=[renderer])

    def _create_legend(self, ax: figure, items: List[Any], title: str) -> None:
        """Create and add a styled legend to the figure."""
        from bokeh.models import Legend

        legend = Legend(
            items=items,
            location="top_right",
            title=title,
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
                self._add_legend_item(ax, source, entry.label, entry.color, marker)
            )
        self._create_legend(ax, items, title or "")

    def hide_spines(self, ax: figure, spines: List[str]) -> None:
        """Hide specified axis spines (no-op for Bokeh).

        Bokeh doesn't have matplotlib-style spines. This method exists
        for interface compatibility but has no visual effect.
        """
        pass

    def hide_yaxis(self, ax: figure) -> None:
        """Hide y-axis ticks, labels, line, and grid for gene track panels."""
        ax.yaxis.visible = False
        ax.ygrid.visible = False

    def format_xaxis_mb(self, ax: figure) -> None:
        """Format x-axis to show megabase values."""
        from bokeh.models import CustomJSTickFormatter

        ax.xaxis.formatter = CustomJSTickFormatter(
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

        Supports .html for interactive and .png/.svg for static.

        Raises:
            ValueError: If the file extension is not supported.
        """
        from bokeh.io import export_png, export_svgs, save
        from bokeh.resources import CDN

        if path.endswith(".html"):
            save(fig, filename=path, resources=CDN)
        elif path.endswith(".png"):
            export_png(fig, filename=path)
        elif path.endswith(".svg"):
            export_svgs(fig, filename=path)
        else:
            raise ValueError(
                f"Unsupported file format: {path!r}. "
                "Supported formats: .html, .png, .svg"
            )

    def show(self, fig: Any) -> None:
        """Display the figure."""
        from bokeh.io import show

        show(fig)

    def close(self, fig: Any) -> None:
        """Close the figure (no-op for bokeh)."""
        pass

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
        from bokeh.models import Whisker

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
        from bokeh.models import BoxAnnotation

        for ax in axes:
            ax.add_layout(
                BoxAnnotation(
                    left=x_start,
                    right=x_end,
                    fill_color=color,
                    fill_alpha=alpha,
                )
            )

    def highlight_heatmap_snp(
        self,
        ax: figure,
        fig: Any,
        snp_idx: int,
        n_snps: int,
        color: str = "#FF0000",
        linewidth: float = 2,
    ) -> None:
        """Highlight a SNP's row and column in the heatmap.

        Uses batched rect() calls for efficiency instead of one renderer
        per cell.

        Args:
            ax: Bokeh figure.
            fig: Layout object (unused, for API compatibility).
            snp_idx: Index of SNP to highlight.
            n_snps: Total number of SNPs in matrix.
            color: Highlight color.
            linewidth: Line width for highlight rectangles.
        """
        if n_snps < 1 or snp_idx < 0 or snp_idx >= n_snps:
            raise ValueError(f"Invalid snp_idx={snp_idx} for n_snps={n_snps}")

        # Collect all highlight cell coordinates, then draw in batched calls
        xs = []
        ys = []

        # Row highlights (lower triangle: columns 0..snp_idx)
        for j in range(snp_idx + 1):
            xs.append(j)
            ys.append(snp_idx)

        # Column highlights (below diagonal: rows snp_idx+1..n_snps-1)
        for i in range(snp_idx + 1, n_snps):
            xs.append(snp_idx)
            ys.append(i)

        if xs:
            source = ColumnDataSource(data={"x": xs, "y": ys})
            ax.rect(
                x="x",
                y="y",
                width=1,
                height=1,
                fill_alpha=0,
                line_color=color,
                line_width=linewidth,
                source=source,
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
        mask_upper: bool = True,
    ) -> Any:
        """Render heatmap with optional triangular masking.

        Args:
            ax: Bokeh figure.
            data: 2D numpy array of values (NaN for missing).
            x_coords: X coordinates for cells.
            y_coords: Y coordinates for cells.
            cmap_colors: Color gradient endpoints [start, end].
            vmin: Minimum value for color scale.
            vmax: Maximum value for color scale.
            mask_upper: If True, mask upper triangle.

        Returns:
            LinearColorMapper for colorbar attachment.
        """
        import numpy as np
        from bokeh.models import LinearColorMapper

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

        # Build rect data (lower triangle only when mask_upper=True)
        xs, ys, values = [], [], []
        n = data.shape[0]
        for i in range(n):
            for j in range(n):
                # Lower triangle including diagonal
                if mask_upper and j > i:
                    continue
                val = data[i, j]
                xs.append(x_coords[j])
                ys.append(y_coords[i])
                values.append(val if not np.isnan(val) else float("nan"))

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

        # Build per-cell width/height arrays matching the flattened xs/ys
        widths, heights = [], []
        for i in range(n):
            for j in range(n):
                if mask_upper and j > i:
                    continue
                widths.append(x_sizes[j])
                heights.append(y_sizes[i])

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
        """Add colorbar legend for heatmap.

        Args:
            ax: Bokeh figure.
            mappable: LinearColorMapper from add_heatmap.
            label: Colorbar label.
            orientation: "vertical" or "horizontal".

        Returns:
            ColorBar object.
        """
        from bokeh.models import BasicTicker, ColorBar

        color_bar = ColorBar(
            color_mapper=mappable,
            ticker=BasicTicker(),
            label_standoff=6,
            title=label,
            orientation=orientation,
        )
        ax.add_layout(color_bar, "right")
        return color_bar


def _parse_color_to_rgb(color: str) -> Tuple[int, int, int]:
    """Parse a color string to an (R, G, B) tuple.

    Supports 6-digit hex (#FF0000), 3-digit hex (#F00), and
    named CSS colors (red, white, etc.) via matplotlib's color converter.

    Args:
        color: Color string in any supported format.

    Returns:
        Tuple of (R, G, B) integers in range 0-255.
    """
    color = color.strip()

    # 6-digit hex
    if color.startswith("#") and len(color) == 7:
        hex_str = color[1:]
        return tuple(int(hex_str[i : i + 2], 16) for i in (0, 2, 4))

    # 3-digit hex — expand to 6-digit
    if color.startswith("#") and len(color) == 4:
        hex_str = color[1:]
        expanded = "".join(c * 2 for c in hex_str)
        return tuple(int(expanded[i : i + 2], 16) for i in (0, 2, 4))

    # Named CSS color — use matplotlib's converter
    try:
        from matplotlib.colors import to_rgb

        r, g, b = to_rgb(color)
        return (int(r * 255), int(g * 255), int(b * 255))
    except (ImportError, ValueError):
        raise ValueError(
            f"Cannot parse color {color!r}. Use 6-digit hex (#FF0000), "
            "3-digit hex (#F00), or a named CSS color (red, blue, etc.)."
        )


def _create_color_palette(start_color: str, end_color: str, n_colors: int) -> List[str]:
    """Create a linear color palette between two colors.

    Args:
        start_color: Starting color (hex or named CSS color).
        end_color: Ending color (hex or named CSS color).
        n_colors: Number of colors in palette.

    Returns:
        List of hex color strings.
    """
    start_rgb = _parse_color_to_rgb(start_color)
    end_rgb = _parse_color_to_rgb(end_color)

    def rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
        return "#{:02x}{:02x}{:02x}".format(*rgb)

    palette = []
    for i in range(n_colors):
        t = i / (n_colors - 1) if n_colors > 1 else 0
        r = int(start_rgb[0] + t * (end_rgb[0] - start_rgb[0]))
        g = int(start_rgb[1] + t * (end_rgb[1] - start_rgb[1]))
        b = int(start_rgb[2] + t * (end_rgb[2] - start_rgb[2]))
        palette.append(rgb_to_hex((r, g, b)))
    return palette

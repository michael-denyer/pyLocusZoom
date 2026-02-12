"""Bokeh backend for pyLocusZoom.

Interactive backend with hover tooltips, well-suited for dashboards.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, DataRange1d, HoverTool, Span
from bokeh.plotting import figure

from . import convert_latex_to_unicode, register_backend

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
    def supports_snp_labels(self) -> bool:
        """Bokeh uses hover tooltips instead of labels."""
        return False

    @property
    def supports_hover(self) -> bool:
        """Bokeh supports hover tooltips."""
        return True

    @property
    def supports_secondary_axis(self) -> bool:
        """Bokeh supports secondary y-axis."""
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

        # Add hover data
        tooltips = []
        if hover_data is not None:
            for col in hover_data.columns:
                data[col] = hover_data[col].values
                col_lower = col.lower()
                if col_lower in ("p-value", "pval", "p_value"):
                    tooltips.append((col, "@{" + col + "}{0.2e}"))
                elif any(x in col_lower for x in ("r2", "r²", "ld")):
                    tooltips.append((col, "@{" + col + "}{0.3f}"))
                elif "pos" in col_lower:
                    tooltips.append((col, "@{" + col + "}{0,0}"))
                else:
                    tooltips.append((col, f"@{col}"))

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
            ax.xaxis.major_label_orientation = (
                rotation * 3.14159 / 180
            )  # Convert to radians

    def _get_legend_location(self, loc: str, default: str = "top_left") -> str:
        """Map matplotlib-style legend location to Bokeh location."""
        loc_map = {
            "upper left": "top_left",
            "upper right": "top_right",
            "lower left": "bottom_left",
            "lower right": "bottom_right",
        }
        return loc_map.get(loc, default)

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

        Returns a dict with configuration for extra_y_ranges.
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

        return "secondary"

    def line_secondary(
        self,
        ax: figure,
        x: pd.Series,
        y: pd.Series,
        color: str = "blue",
        linewidth: float = 1.5,
        alpha: float = 1.0,
        linestyle: str = "-",
        label: Optional[str] = None,
        yaxis_name: str = "secondary",
    ) -> Any:
        """Create a line plot on secondary y-axis."""
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
        ax: figure,
        x: pd.Series,
        y1: Union[float, pd.Series],
        y2: Union[float, pd.Series],
        color: str = "blue",
        alpha: float = 0.3,
        yaxis_name: str = "secondary",
    ) -> Any:
        """Fill area between two y-values on secondary y-axis."""
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
        ax: figure,
        bottom: float,
        top: float,
        yaxis_name: str = "secondary",
    ) -> None:
        """Set secondary y-axis limits."""
        if yaxis_name in ax.extra_y_ranges:
            ax.extra_y_ranges[yaxis_name].start = bottom
            ax.extra_y_ranges[yaxis_name].end = top

    def set_secondary_ylabel(
        self,
        ax: figure,
        label: str,
        color: str = "black",
        fontsize: int = 10,
        yaxis_name: str = "secondary",
    ) -> None:
        """Set secondary y-axis label."""
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

    def add_snp_labels(
        self,
        ax: figure,
        df: pd.DataFrame,
        pos_col: str,
        neglog10p_col: str,
        rs_col: str,
        label_top_n: int,
        genes_df: Optional[pd.DataFrame],
        chrom: int,
        adjust: bool = True,
    ) -> List[Any]:
        """No-op: Bokeh uses hover tooltips instead of text labels."""
        return []

    def adjust_snp_labels(self, ax: Any, texts: List[Any]) -> None:
        """No-op: Bokeh uses hover tooltips instead of text labels."""
        pass

    def add_panel_label(
        self,
        ax: figure,
        label: str,
        x_frac: float = 0.02,
        y_frac: float = 0.95,
    ) -> None:
        """Add label text at fractional position in panel."""
        from bokeh.models import Label

        # Convert fraction to data coordinates using axis ranges
        x_range = ax.x_range
        y_range = ax.y_range
        x = (
            x_range.start + x_frac * (x_range.end - x_range.start)
            if hasattr(x_range, "start") and x_range.start is not None
            else 0
        )

        # Handle both normal and inverted y-axis
        # For inverted axis (start > end), y_frac=0.95 should be near start (top visually)
        if hasattr(y_range, "start") and y_range.start is not None:
            y_min = min(y_range.start, y_range.end)
            y_max = max(y_range.start, y_range.end)
            # Position label at top of visible range regardless of axis direction
            y = y_min + y_frac * (y_max - y_min)
        else:
            y = 0

        label_obj = Label(
            x=x,
            y=y,
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
        from bokeh.models import ColumnDataSource, Range1d

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

    def add_ld_legend(
        self,
        ax: figure,
        ld_bins: List[Tuple[float, str, str]],
        lead_snp_color: str,
    ) -> None:
        """Add LD color legend using invisible dummy glyphs.

        Creates legend entries with dummy renderers that are excluded from
        the data range calculation to avoid affecting axis scaling.
        """
        source = self._ensure_legend_range(ax)
        items = [
            self._add_legend_item(ax, source, "Lead SNP", lead_snp_color, "diamond", 16)
        ]
        for _, label, color in ld_bins:
            items.append(self._add_legend_item(ax, source, label, color, "square"))
        self._create_legend(ax, items, "r²")

    def add_effect_legend(
        self,
        ax: figure,
        effect_bins: List[Tuple[float, str, str]],
    ) -> None:
        """Add effect direction legend for colocalization plots."""
        source = self._ensure_legend_range(ax)
        items = []
        for _, label, color in effect_bins:
            items.append(self._add_legend_item(ax, source, label, color, "circle"))
        self._create_legend(ax, items, "Effect")

    def add_legend(
        self,
        ax: figure,
        handles: List[Any],
        labels: List[str],
        loc: str = "upper left",
        title: Optional[str] = None,
    ) -> Any:
        """Configure legend on the figure."""
        ax.legend.location = self._get_legend_location(loc, "top_left")
        if title:
            ax.legend.title = title
        ax.legend.background_fill_alpha = 0.9
        ax.legend.border_line_color = "black"

        return ax.legend

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

        Supports .html for interactive and .png for static.
        """
        from bokeh.io import export_png, export_svgs, output_file, save

        if path.endswith(".html"):
            output_file(path)
            save(fig)
        elif path.endswith(".png"):
            export_png(fig, filename=path)
        elif path.endswith(".svg"):
            export_svgs(fig, filename=path)
        else:
            output_file(path)
            save(fig)

    def show(self, fig: Any) -> None:
        """Display the figure."""
        from bokeh.io import show

        show(fig)

    def close(self, fig: Any) -> None:
        """Close the figure (no-op for bokeh)."""
        pass

    def add_eqtl_legend(
        self,
        ax: figure,
        eqtl_positive_bins: List[Tuple[float, float, str, str]],
        eqtl_negative_bins: List[Tuple[float, float, str, str]],
    ) -> None:
        """Add eQTL effect size legend using invisible dummy glyphs."""
        source = self._ensure_legend_range(ax)
        items = []
        for _, _, label, color in eqtl_positive_bins:
            items.append(self._add_legend_item(ax, source, label, color, "triangle"))
        for _, _, label, color in eqtl_negative_bins:
            items.append(
                self._add_legend_item(ax, source, label, color, "inverted_triangle")
            )
        self._create_legend(ax, items, "eQTL effect")

    def add_finemapping_legend(
        self,
        ax: figure,
        credible_sets: List[int],
        get_color_func: Any,
    ) -> None:
        """Add fine-mapping credible set legend using invisible dummy glyphs."""
        if not credible_sets:
            return

        source = self._ensure_legend_range(ax)
        items = [
            self._add_legend_item(
                ax, source, f"CS{cs_id}", get_color_func(cs_id), "circle"
            )
            for cs_id in credible_sets
        ]
        self._create_legend(ax, items, "Credible sets")

    def add_simple_legend(
        self,
        ax: figure,
        label: str,
        loc: str = "upper right",
    ) -> None:
        """Configure legend position.

        Bokeh handles legends automatically from legend_label.
        This just positions the legend.
        """
        ax.legend.location = self._get_legend_location(loc, "top_right")
        ax.legend.background_fill_alpha = 0.9
        ax.legend.border_line_color = "black"

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

    def add_recombination_overlay(
        self,
        ax: figure,
        recomb_df: pd.DataFrame,
        start: int,
        end: int,
    ) -> None:
        """Add recombination overlay using bokeh secondary y-axis.

        Args:
            ax: Bokeh figure.
            recomb_df: DataFrame with 'pos' and 'rate' columns.
            start: Region start position for filtering.
            end: Region end position for filtering.
        """
        from ..recombination import RECOMB_COLOR

        # Filter to region
        region_recomb = recomb_df[
            (recomb_df["pos"] >= start) & (recomb_df["pos"] <= end)
        ].copy()

        if region_recomb.empty:
            return

        # Create twin axis - returns "secondary" string
        yaxis_name = self.create_twin_axis(ax)

        # Plot fill under curve
        self.fill_between_secondary(
            ax,
            region_recomb["pos"],
            0,
            region_recomb["rate"],
            color=RECOMB_COLOR,
            alpha=0.15,
            yaxis_name=yaxis_name,
        )

        # Plot recombination rate line
        self.line_secondary(
            ax,
            region_recomb["pos"],
            region_recomb["rate"],
            color=RECOMB_COLOR,
            linewidth=2.5,
            alpha=0.8,
            yaxis_name=yaxis_name,
        )

        # Set y-axis limits and label
        max_rate = region_recomb["rate"].max()
        self.set_secondary_ylim(ax, 0, max(max_rate * 1.3, 10), yaxis_name=yaxis_name)
        self.set_secondary_ylabel(
            ax,
            "Recombination rate (cM/Mb)",
            color="black",
            fontsize=9,
            yaxis_name=yaxis_name,
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

        # Build rect data for lower triangle
        xs, ys, values = [], [], []
        n = data.shape[0]
        for i in range(n):
            for j in range(n):
                # Lower triangle including diagonal
                if mask_upper and j > i:
                    continue
                val = data[i, j]
                xs.append(j)
                ys.append(i)
                values.append(val if not np.isnan(val) else float("nan"))

        source = ColumnDataSource({"x": xs, "y": ys, "value": values})

        ax.rect(
            x="x",
            y="y",
            width=1,
            height=1,
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

        orientation_map = {"vertical": "vertical", "horizontal": "horizontal"}
        color_bar = ColorBar(
            color_mapper=mappable,
            ticker=BasicTicker(),
            label_standoff=6,
            title=label,
            orientation=orientation_map.get(orientation, "vertical"),
        )
        ax.add_layout(color_bar, "right")
        return color_bar


def _create_color_palette(start_color: str, end_color: str, n_colors: int) -> List[str]:
    """Create a linear color palette between two colors.

    Args:
        start_color: Starting hex color.
        end_color: Ending hex color.
        n_colors: Number of colors in palette.

    Returns:
        List of hex color strings.
    """

    def hex_to_rgb(hex_color: str) -> Tuple[int, int, int]:
        hex_color = hex_color.lstrip("#")
        return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    def rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
        return "#{:02x}{:02x}{:02x}".format(*rgb)

    start_rgb = hex_to_rgb(start_color)
    end_rgb = hex_to_rgb(end_color)

    palette = []
    for i in range(n_colors):
        t = i / (n_colors - 1) if n_colors > 1 else 0
        r = int(start_rgb[0] + t * (end_rgb[0] - start_rgb[0]))
        g = int(start_rgb[1] + t * (end_rgb[1] - start_rgb[1]))
        b = int(start_rgb[2] + t * (end_rgb[2] - start_rgb[2]))
        palette.append(rgb_to_hex((r, g, b)))
    return palette

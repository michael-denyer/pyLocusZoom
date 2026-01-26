"""Bokeh backend for pyLocusZoom.

Interactive backend with hover tooltips, well-suited for dashboards.
"""

from typing import Any, List, Optional, Tuple, Union

import pandas as pd
from bokeh.io import export_png, export_svgs, output_file, save, show
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, DataRange1d, HoverTool, Span
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
        x_range = DataRange1d() if sharex else None

        for i, h in enumerate(heights):
            p = figure(
                width=width_px,
                height=h,
                x_range=x_range if sharex else DataRange1d(),
                tools="pan,wheel_zoom,box_zoom,reset,save",
                toolbar_location="above" if i == 0 else None,
            )

            # Style
            p.grid.grid_line_alpha = 0.3
            p.outline_line_color = None

            figures.append(p)

        # Create column layout (use default sizing mode to avoid validation warnings)
        layout = column(*figures)

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

        # Get marker type for scatter()
        marker_type = self._marker_map.get(marker, "circle")

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
        # Convert linestyle
        dash_map = {
            "-": "solid",
            "--": "dashed",
            ":": "dotted",
            "-.": "dashdot",
        }
        line_dash = dash_map.get(linestyle, "solid")

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
        dash_map = {"-": "solid", "--": "dashed", ":": "dotted", "-.": "dashdot"}
        line_dash = dash_map.get(linestyle, "dashed")

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
        label = label.replace("$", "")
        return label

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
        dash_map = {"-": "solid", "--": "dashed", ":": "dotted", "-.": "dashdot"}
        line_dash = dash_map.get(linestyle, "solid")

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
        from bokeh.models import ColumnDataSource, Legend, LegendItem, Range1d, Scatter

        legend_items = []

        # Create a separate range for legend glyphs that won't affect the main plot
        if "legend_range" not in ax.extra_y_ranges:
            ax.extra_y_ranges["legend_range"] = Range1d(start=0, end=1)

        # Use coordinates within the legend range
        dummy_source = ColumnDataSource(data={"x": [0], "y": [0]})

        # Add LD bin markers (no lead SNP - it's shown in the actual plot)
        for _, label, color in ld_bins:
            glyph = Scatter(
                x="x",
                y="y",
                marker="square",
                size=10,
                fill_color=color,
                line_color="black",
                line_width=0.5,
            )
            renderer = ax.add_glyph(dummy_source, glyph)
            renderer.y_range_name = "legend_range"
            renderer.visible = False
            legend_items.append(LegendItem(label=label, renderers=[renderer]))

        legend = Legend(
            items=legend_items,
            location="top_right",
            title="r²",
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

    def add_eqtl_legend(
        self,
        ax: figure,
        eqtl_positive_bins: List[Tuple[float, float, str, str]],
        eqtl_negative_bins: List[Tuple[float, float, str, str]],
    ) -> None:
        """Add eQTL effect size legend using invisible dummy glyphs."""
        from bokeh.models import ColumnDataSource, Legend, LegendItem, Range1d

        legend_items = []

        # Create a separate range for legend glyphs that won't affect the main plot
        if "legend_range" not in ax.extra_y_ranges:
            ax.extra_y_ranges["legend_range"] = Range1d(start=0, end=1)

        dummy_source = ColumnDataSource(data={"x": [0], "y": [0]})

        # Positive effects (upward triangles)
        for _, _, label, color in eqtl_positive_bins:
            renderer = ax.scatter(
                x="x",
                y="y",
                source=dummy_source,
                marker="triangle",
                size=10,
                fill_color=color,
                line_color="black",
                line_width=0.5,
                y_range_name="legend_range",
                visible=False,
            )
            legend_items.append(LegendItem(label=label, renderers=[renderer]))

        # Negative effects (downward triangles)
        for _, _, label, color in eqtl_negative_bins:
            renderer = ax.scatter(
                x="x",
                y="y",
                source=dummy_source,
                marker="inverted_triangle",
                size=10,
                fill_color=color,
                line_color="black",
                line_width=0.5,
                y_range_name="legend_range",
                visible=False,
            )
            legend_items.append(LegendItem(label=label, renderers=[renderer]))

        legend = Legend(
            items=legend_items,
            location="top_right",
            title="eQTL effect",
            background_fill_alpha=0.9,
            border_line_color="black",
            spacing=0,
            padding=4,
            label_height=12,
            glyph_height=12,
        )
        ax.add_layout(legend)

    def add_finemapping_legend(
        self,
        ax: figure,
        credible_sets: List[int],
        get_color_func: Any,
    ) -> None:
        """Add fine-mapping credible set legend using invisible dummy glyphs."""
        if not credible_sets:
            return

        from bokeh.models import ColumnDataSource, Legend, LegendItem, Range1d

        legend_items = []

        # Create a separate range for legend glyphs that won't affect the main plot
        if "legend_range" not in ax.extra_y_ranges:
            ax.extra_y_ranges["legend_range"] = Range1d(start=0, end=1)

        dummy_source = ColumnDataSource(data={"x": [0], "y": [0]})

        for cs_id in credible_sets:
            color = get_color_func(cs_id)
            renderer = ax.scatter(
                x="x",
                y="y",
                source=dummy_source,
                marker="circle",
                size=10,
                fill_color=color,
                line_color="black",
                line_width=0.5,
                y_range_name="legend_range",
                visible=False,
            )
            legend_items.append(LegendItem(label=f"CS{cs_id}", renderers=[renderer]))

        legend = Legend(
            items=legend_items,
            location="top_right",
            title="Credible sets",
            background_fill_alpha=0.9,
            border_line_color="black",
            spacing=0,
            padding=4,
            label_height=12,
            glyph_height=12,
        )
        ax.add_layout(legend)

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

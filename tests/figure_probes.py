"""Per-backend inspection of a rendered figure, behind one shared vocabulary.

Matplotlib, Plotly and Bokeh answer the same questions in different words:
"what marker symbols are on this figure" is ``trace.marker.symbol`` in plotly
and a walk over ``figure.renderers[*].glyph.marker`` in bokeh, and each library
names the symbols differently. A probe translates its library's answer into the
shared vocabulary below, so a test can state the behaviour once and run it on
every backend.

These three classes are the only place in the suite that knows the libraries'
figure internals. A behaviour test that reaches past a probe is a test that
will drift from its twin.

Questions that only an interactive figure can answer (``marker_symbols``,
``has_hover``, ``hover_values``, ``standalone_html``, ``json_payload``) are
asked over ``INTERACTIVE_BACKENDS``, so ``MatplotlibProbe`` does not define them.

Vocabulary:

- Markers, as plotly already spells them: ``circle``, ``diamond``,
  ``triangle-up``, ``triangle-down``.
- Colours: lower-case ``#rrggbb`` hex.
- Legend corners, as matplotlib spells ``loc``: ``"lower left"``.
- Panels: the stacked plotting areas, top first, as ``create_figure`` returned
  them. A matplotlib twin axis or colorbar is not a panel.
- ``Box``: an axis-aligned rectangle drawn in data coordinates with
  ``add_rectangle``. A ``RegionHighlight`` spans the panel's full height.
"""

import json
from typing import NamedTuple, Optional

from matplotlib.colors import to_hex

BOKEH_MARKER_NAMES = {
    "triangle": "triangle-up",
    "inverted_triangle": "triangle-down",
}

INTERACTIVE_BACKENDS = ("plotly", "bokeh")

LEGEND_VERTICAL = {"bottom": "lower", "top": "upper"}


class Box(NamedTuple):
    """A rectangle drawn in data coordinates."""

    x0: float
    x1: float
    y0: float
    y1: float
    facecolor: Optional[str]
    edgecolor: Optional[str]

    @property
    def centre(self):
        return ((self.x0 + self.x1) / 2, (self.y0 + self.y1) / 2)


class RegionHighlight(NamedTuple):
    """A full-height band shading an x-range of a panel."""

    x0: float
    x1: float
    color: str


def _hex(color):
    """Normalise any colour spelling to lower-case hex, None for no colour."""
    if color is None:
        return None
    if isinstance(color, str) and color.replace(" ", "") == "rgba(0,0,0,0)":
        return None
    return to_hex(color, keep_alpha=False)


def _pt(size):
    """Font size in points from a number or a ``"12pt"`` string."""
    return float(str(size).removesuffix("pt"))


class MatplotlibProbe:
    """Read a matplotlib Figure."""

    @staticmethod
    def panels(fig):
        """The stacked panels, dropping twin axes and colorbars."""
        seen, panels = set(), []
        for ax in fig.get_axes():
            if ax.get_subplotspec() is None:
                continue
            bounds = tuple(round(b, 6) for b in ax.get_position().bounds)
            if bounds not in seen:
                seen.add(bounds)
                panels.append(ax)
        return panels

    def panel_count(self, fig):
        """How many stacked panels the figure carries."""
        return len(self.panels(fig))

    def xticks(self, fig, panel=0):
        """Tick positions and labels on one panel's x-axis."""
        ax = self.panels(fig)[panel]
        ticks = [float(t) for t in ax.get_xticks()]
        formatter = ax.xaxis.get_major_formatter()
        return ticks, [formatter(t, i) for i, t in enumerate(ticks)]

    def x_axis_in_mb(self, fig, panel=0):
        """Whether the panel labels x positions in megabases."""
        formatter = self.panels(fig)[panel].xaxis.get_major_formatter()
        return formatter(1_500_000, 0) == "1.50"

    def legend_corner(self, fig, panel=0):
        """The corner the panel's legend is anchored to."""
        from matplotlib.legend import Legend

        legend = self.panels(fig)[panel].get_legend()
        names = {code: name for name, code in Legend.codes.items()}
        return names[legend._get_loc()]

    def legend_edgecolors(self, fig, panel=0):
        """Each legend label mapped to its swatch's edge colour."""
        legend = self.panels(fig)[panel].get_legend()
        edges = {}
        for text, handle in zip(legend.get_texts(), legend.legend_handles):
            edge = (
                handle.get_markeredgecolor()
                if hasattr(handle, "get_markeredgecolor")
                else handle.get_edgecolor()
            )
            edges[text.get_text()] = _hex(edge)
        return edges

    def hline_levels(self, fig, panel=None, linestyle="--"):
        """Heights of full-width horizontal lines, in panel order."""
        panels = self.panels(fig) if panel is None else [self.panels(fig)[panel]]
        return [
            float(line.get_ydata()[0])
            for ax in panels
            for line in ax.get_lines()
            if list(line.get_xdata()) == [0, 1]
            and (linestyle is None or line.get_linestyle() == linestyle)
        ]

    def marker_x(self, fig, panel=0, color=None):
        """Sorted x of every scatter marker on one panel, optionally of one fill."""
        from matplotlib.collections import PathCollection

        xs = []
        for c in self.panels(fig)[panel].collections:
            if not isinstance(c, PathCollection):
                continue
            offsets = c.get_offsets()
            fills = [_hex(f) for f in c.get_facecolors()] or [None]
            if len(fills) == 1:
                fills = fills * len(offsets)
            xs.extend(
                float(x)
                for (x, _), fill in zip(offsets, fills)
                if color is None or fill == _hex(color)
            )
        return sorted(xs)

    def point_alphas(self, fig):
        """The alpha of every scatter layer."""
        from matplotlib.collections import PathCollection

        return {
            c.get_alpha()
            for ax in self.panels(fig)
            for c in ax.collections
            if isinstance(c, PathCollection)
        }

    def boxes(self, fig, panel=0):
        """Rectangles drawn in data coordinates on one panel."""
        from matplotlib.patches import Rectangle

        ax = self.panels(fig)[panel]
        return [
            Box(
                p.get_x(),
                p.get_x() + p.get_width(),
                p.get_y(),
                p.get_y() + p.get_height(),
                _hex(p.get_facecolor()) if p.get_fill() else None,
                _hex(p.get_edgecolor()),
            )
            for p in ax.patches
            if isinstance(p, Rectangle) and p.get_data_transform() is ax.transData
        ]

    def region_highlights(self, fig, panel=0):
        """Full-height bands shading an x-range of one panel."""
        from matplotlib.patches import Rectangle

        ax = self.panels(fig)[panel]
        return [
            RegionHighlight(
                p.get_x(), p.get_x() + p.get_width(), _hex(p.get_facecolor())
            )
            for p in ax.patches
            if isinstance(p, Rectangle) and p.get_data_transform() is not ax.transData
        ]

    def colorbar_titles(self, fig):
        """The title of every colour scale on the figure."""
        return [ax.get_ylabel() for ax in fig.axes if hasattr(ax, "_colorbar")]

    def font_sizes(self, fig):
        """Point sizes of the figure title, panel titles, axis labels and ticks."""
        panels = self.panels(fig)
        return {
            "title": fig._suptitle.get_fontsize() if fig._suptitle else None,
            "panel_title": {ax.title.get_fontsize() for ax in panels if ax.get_title()},
            "axis_label": {
                label.get_fontsize()
                for ax in panels
                for label in (ax.xaxis.label, ax.yaxis.label)
            },
            "tick_label": {
                t.get_fontsize()
                for ax in panels
                for t in ax.get_xticklabels() + ax.get_yticklabels()
            },
        }


class PlotlyProbe:
    """Read a plotly Figure."""

    @staticmethod
    def _axis_keys(fig, kind):
        layout = fig.layout.to_plotly_json()
        return sorted(
            (
                key
                for key in layout
                if key.startswith(kind) and "overlaying" not in layout[key]
            ),
            key=lambda k: int(k.removeprefix(kind) or 1),
        )

    def panel_count(self, fig):
        """How many subplots the figure carries."""
        return len(self._axis_keys(fig, "xaxis"))

    def _xref(self, fig, panel):
        """The x reference a shape on ``panel`` carries, such as ``"x2"``."""
        key = self._axis_keys(fig, "xaxis")[panel]
        return "x" + key.removeprefix("xaxis")

    def _xaxis(self, fig, panel):
        return fig.layout[self._axis_keys(fig, "xaxis")[panel]]

    def marker_symbols(self, fig):
        """Every marker symbol drawn, in the shared vocabulary."""
        return {
            str(trace.marker.symbol)
            for trace in fig.data
            if getattr(trace, "marker", None) is not None
            and trace.marker.symbol is not None
        }

    def has_hover(self, fig):
        """Whether any trace carries hover text."""
        return any(getattr(trace, "hovertemplate", None) for trace in fig.data)

    def hover_values(self, fig):
        """Every value a hover tooltip can show, as strings."""
        return {
            str(value)
            for trace in fig.data
            if getattr(trace, "customdata", None) is not None
            for row in trace.customdata
            for value in row
        }

    def standalone_html(self, fig):
        """A complete HTML document carrying the figure and its library."""
        return fig.to_html(include_plotlyjs=True, full_html=True)

    def json_payload(self, fig):
        """The figure as the dict a notebook front end is handed."""
        return json.loads(fig.to_json())

    def heatmap_coords(self, ax):
        """The x and y coordinates of the heatmap drawn on one panel."""
        trace = next(t for t in ax[0].data if hasattr(t, "z"))
        return list(trace.x), list(trace.y)

    def xticks(self, fig, panel=0):
        """Tick positions and labels on one panel's x-axis."""
        axis = self._xaxis(fig, panel)
        return list(axis.tickvals or []), list(axis.ticktext or [])

    def x_axis_in_mb(self, fig, panel=0):
        """Whether the panel labels x positions in megabases."""
        axis = self._xaxis(fig, panel)
        if axis.ticksuffix != " Mb" or not axis.tickvals:
            return False
        return [f"{v / 1e6:.2f}" for v in axis.tickvals] == list(axis.ticktext)

    def legend_corner(self, fig, panel=0):
        """The corner the first legend is anchored to."""
        legend = fig.layout.legend
        return f"{LEGEND_VERTICAL[legend.yanchor]} {legend.xanchor}"

    def legend_edgecolors(self, fig, panel=0):
        """Each legend label mapped to its swatch's edge colour."""
        return {
            trace.name: _hex(trace.marker.line.color)
            for trace in fig.data
            if trace.showlegend and trace.name
        }

    def _shapes(self, fig, panel, shape_type):
        xref = None if panel is None else self._xref(fig, panel)
        return [
            s
            for s in fig.layout.shapes
            if s.type == shape_type
            and (xref is None or s.xref in (xref, f"{xref} domain"))
        ]

    def hline_levels(self, fig, panel=None, linestyle="--"):
        """Heights of full-width horizontal lines, in panel order."""
        dash = {"--": "dash", ":": "dot", "-": "solid", None: None}[linestyle]
        return [
            float(s.y0)
            for s in self._shapes(fig, panel, "line")
            if str(s.xref).endswith("domain")
            and s.y0 == s.y1
            and (dash is None or s.line.dash == dash)
        ]

    def marker_x(self, fig, panel=0, color=None):
        """Sorted x of every scatter marker on one panel, optionally of one fill."""
        xref = self._xref(fig, panel)
        xs = []
        for t in fig.data:
            if t.type != "scatter" or t.mode != "markers" or t.x is None:
                continue
            if (t.xaxis or "x") != xref:
                continue
            fills = t.marker.color
            if fills is None or isinstance(fills, str):
                fills = [fills] * len(t.x)
            xs.extend(
                float(x)
                for x, fill in zip(t.x, fills)
                if x is not None and (color is None or _hex(fill) == _hex(color))
            )
        return sorted(xs)

    def point_alphas(self, fig):
        """The opacity of every marker trace that carries data."""
        return {
            t.marker.opacity
            for t in fig.data
            if t.type == "scatter"
            and t.mode == "markers"
            and t.x is not None
            and any(x is not None for x in list(t.x))
        }

    def boxes(self, fig, panel=0):
        """Rectangles drawn in data coordinates on one panel."""
        return [
            Box(s.x0, s.x1, s.y0, s.y1, _hex(s.fillcolor), _hex(s.line.color))
            for s in self._shapes(fig, panel, "rect")
            if not str(s.yref).endswith("domain")
        ]

    def region_highlights(self, fig, panel=0):
        """Full-height bands shading an x-range of one panel."""
        return [
            RegionHighlight(s.x0, s.x1, _hex(s.fillcolor))
            for s in self._shapes(fig, panel, "rect")
            if str(s.yref).endswith("domain")
        ]

    def colorbar_titles(self, fig):
        """The title of every colour scale on the figure."""
        return [
            trace.colorbar.title.text
            for trace in fig.data
            if trace.type == "heatmap" and trace.showscale
        ]

    def font_sizes(self, fig):
        """Point sizes of the figure title, panel titles, axis labels and ticks."""
        layout = fig.layout
        axes = [
            layout[k] for kind in ("xaxis", "yaxis") for k in self._axis_keys(fig, kind)
        ]
        return {
            "title": layout.title.font.size,
            "panel_title": {a.font.size for a in layout.annotations},
            "axis_label": {a.title.font.size for a in axes},
            "tick_label": {a.tickfont.size for a in axes},
        }


class BokehProbe:
    """Read a bokeh LayoutDOM."""

    @staticmethod
    def panels(fig):
        """The plots in the layout, top to bottom, left to right."""
        from bokeh.models import Plot

        if isinstance(fig, Plot):
            return [fig]
        return [
            p
            for child in getattr(fig, "children", [])
            for p in BokehProbe.panels(child)
        ]

    def panel_count(self, fig):
        """How many plots the layout carries."""
        return len(self.panels(fig))

    def _scatter_glyphs(self, fig):
        from bokeh.models import GlyphRenderer, Scatter

        return [
            renderer.glyph
            for plot in self.panels(fig)
            for renderer in plot.renderers
            if isinstance(renderer, GlyphRenderer)
            and isinstance(renderer.glyph, Scatter)
        ]

    def marker_symbols(self, fig):
        """Every marker symbol drawn, in the shared vocabulary."""
        return {
            BOKEH_MARKER_NAMES.get(glyph.marker, glyph.marker)
            for glyph in self._scatter_glyphs(fig)
        }

    def has_hover(self, fig):
        """Whether any panel carries a hover tool."""
        from bokeh.models import HoverTool

        return any(
            isinstance(tool, HoverTool)
            for plot in self.panels(fig)
            for tool in plot.tools
        )

    def hover_values(self, fig):
        """Every value a hover tooltip can show, as strings."""
        import re

        from bokeh.models import GlyphRenderer, HoverTool

        values = set()
        for plot in self.panels(fig):
            fields = {
                braced or bare
                for tool in plot.tools
                if isinstance(tool, HoverTool)
                for _, spec in tool.tooltips or []
                for braced, bare in re.findall(r"@(?:\{([^}]+)\}|(\w+))", spec)
            }
            for renderer in plot.renderers:
                if isinstance(renderer, GlyphRenderer):
                    data = renderer.data_source.data
                    for field in fields & set(data):
                        values.update(str(v) for v in data[field])
        return values

    def standalone_html(self, fig):
        """A complete HTML document carrying the figure and its library."""
        from bokeh.embed import file_html
        from bokeh.resources import CDN

        return file_html(fig, CDN)

    def json_payload(self, fig):
        """The figure as the dict a notebook front end is handed."""
        from bokeh.embed import json_item

        return json_item(fig)

    def heatmap_coords(self, ax):
        """The x and y coordinates of the heatmap drawn on one panel."""
        rect = [r for r in ax.renderers if hasattr(r, "glyph")][-1]
        source = rect.data_source.data
        return sorted(set(source["x"])), sorted(set(source["y"]))

    def xticks(self, fig, panel=0):
        """Tick positions and labels on one panel's x-axis."""
        axis = self.panels(fig)[panel].xaxis[0]
        ticks = list(axis.ticker.ticks)
        overrides = axis.major_label_overrides
        return ticks, [str(overrides.get(t, t)) for t in ticks] if overrides else []

    def x_axis_in_mb(self, fig, panel=0):
        """Whether the panel labels x positions in megabases."""
        from bokeh.models import CustomJSTickFormatter

        formatter = self.panels(fig)[panel].xaxis[0].formatter
        return isinstance(formatter, CustomJSTickFormatter) and "1e6" in formatter.code

    def _legend(self, fig, panel):
        from bokeh.models import Legend

        return next(iter(self.panels(fig)[panel].select(Legend)))

    def legend_corner(self, fig, panel=0):
        """The corner the panel's legend is anchored to."""
        vertical, horizontal = self._legend(fig, panel).location.split("_")
        return f"{LEGEND_VERTICAL[vertical]} {horizontal}"

    def legend_edgecolors(self, fig, panel=0):
        """Each legend label mapped to its swatch's edge colour."""
        return {
            item.label.value: _hex(item.renderers[0].glyph.line_color)
            for item in self._legend(fig, panel).items
        }

    def hline_levels(self, fig, panel=None, linestyle="--"):
        """Heights of full-width horizontal lines, in panel order."""
        from bokeh.models import Span

        dash = {"--": "dashed", ":": "dotted", "-": "solid", None: None}[linestyle]
        plots = self.panels(fig) if panel is None else [self.panels(fig)[panel]]
        return [
            float(span.location)
            for plot in plots
            for span in plot.center
            if isinstance(span, Span)
            and span.dimension == "width"
            and (dash is None or _dash_name(span.line_dash) == dash)
        ]

    def marker_x(self, fig, panel=0, color=None):
        """Sorted x of every scatter marker on one panel, optionally of one fill."""
        from bokeh.models import GlyphRenderer, Scatter

        xs = []
        for renderer in self.panels(fig)[panel].renderers:
            if not (
                isinstance(renderer, GlyphRenderer)
                and isinstance(renderer.glyph, Scatter)
            ):
                continue
            positions = _glyph_values(renderer, "x")
            fills = _glyph_values(renderer, "fill_color")
            xs.extend(
                float(x)
                for x, fill in zip(positions, fills)
                if color is None or _hex(fill) == _hex(color)
            )
        return sorted(xs)

    def point_alphas(self, fig):
        """The alpha of every scatter glyph, where fill and line agree."""
        return {
            glyph.fill_alpha if glyph.fill_alpha == glyph.line_alpha else None
            for glyph in self._scatter_glyphs(fig)
        }

    def boxes(self, fig, panel=0):
        """Rectangles drawn in data coordinates on one panel."""
        from bokeh.models import GlyphRenderer, Rect

        boxes = []
        for renderer in self.panels(fig)[panel].renderers:
            if not isinstance(renderer, GlyphRenderer):
                continue
            glyph = renderer.glyph
            if not isinstance(glyph, Rect) or "value" in renderer.data_source.data:
                continue
            columns = [
                _glyph_values(renderer, p) for p in ("x", "y", "width", "height")
            ]
            for x, y, w, h in zip(*columns):
                boxes.append(
                    Box(
                        x - w / 2,
                        x + w / 2,
                        y - h / 2,
                        y + h / 2,
                        _hex(glyph.fill_color),
                        _hex(glyph.line_color),
                    )
                )
        return boxes

    def region_highlights(self, fig, panel=0):
        """Full-height bands shading an x-range of one panel."""
        from bokeh.models import BoxAnnotation

        return [
            RegionHighlight(box.left, box.right, _hex(box.fill_color))
            for box in self.panels(fig)[panel].center
            if isinstance(box, BoxAnnotation)
        ]

    def colorbar_titles(self, fig):
        """The title of every colour scale on the figure."""
        from bokeh.models import ColorBar

        return [
            bar.title
            for plot in self.panels(fig)
            for bar in plot.right
            if isinstance(bar, ColorBar)
        ]

    def font_sizes(self, fig):
        """Point sizes of the figure title, panel titles, axis labels and ticks."""
        from bokeh.models import Div

        plots = self.panels(fig)
        title = next(
            (c for c in getattr(fig, "children", []) if isinstance(c, Div)), None
        )
        axes = [axis for plot in plots for axis in plot.xaxis + plot.yaxis]
        return {
            "title": _pt(title.styles["font-size"]) if title else None,
            "panel_title": {_pt(p.title.text_font_size) for p in plots if p.title.text},
            "axis_label": {_pt(a.axis_label_text_font_size) for a in axes},
            "tick_label": {_pt(a.major_label_text_font_size) for a in axes},
        }


def _glyph_values(renderer, prop):
    """Per-item values of a bokeh glyph property, whether column or literal."""
    spec = getattr(renderer.glyph, prop)
    data = renderer.data_source.data
    if isinstance(spec, str) and spec in data:
        return list(data[spec])
    if isinstance(spec, dict) and "field" in spec:
        return list(data[spec["field"]])
    return [spec] * len(next(iter(data.values())))


def _dash_name(line_dash):
    """Bokeh stores a named dash as its pattern; name the common ones."""
    patterns = {(): "solid", (6,): "dashed", (2, 4): "dotted"}
    if isinstance(line_dash, str):
        return line_dash
    return patterns.get(tuple(line_dash), "custom")


PROBES = {
    "matplotlib": MatplotlibProbe(),
    "plotly": PlotlyProbe(),
    "bokeh": BokehProbe(),
}

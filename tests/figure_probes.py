"""Per-backend inspection of a rendered figure, behind one shared vocabulary.

Plotly and Bokeh answer the same questions in different words: "what marker
symbols are on this figure" is ``trace.marker.symbol`` in one and a walk over
``figure.renderers[*].glyph.marker`` in the other, and each library names the
symbols differently. A probe translates its library's answer into the shared
vocabulary below, so a test can state the behaviour once and run it on every
interactive backend.

These two classes are the only place in the suite that knows either library's
internals. A behaviour test that reaches past a probe is a test that will
drift from its twin.

Matplotlib is absent on purpose. It has no hover layer and no self-contained
HTML export, so the questions here have no answer for it; the matplotlib
figures are inspected directly by ``test_plotter_backends.py``.

Marker vocabulary, as plotly already spells it: ``circle``, ``diamond``,
``triangle-up``, ``triangle-down``.
"""

import json

BOKEH_MARKER_NAMES = {
    "triangle": "triangle-up",
    "inverted_triangle": "triangle-down",
}

INTERACTIVE_BACKENDS = ("plotly", "bokeh")


class PlotlyProbe:
    """Read a plotly Figure."""

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


class BokehProbe:
    """Read a bokeh LayoutDOM."""

    @staticmethod
    def _scatter_glyphs(fig):
        from bokeh.models import GlyphRenderer, Scatter

        return [
            renderer.glyph
            for child in fig.children
            for renderer in getattr(child, "renderers", [])
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
            for child in fig.children
            for tool in getattr(child, "tools", [])
        )

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


PROBES = {"plotly": PlotlyProbe(), "bokeh": BokehProbe()}

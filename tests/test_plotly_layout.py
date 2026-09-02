"""Tests for the plotly subplot arithmetic in backends/plotly_layout.py.

These functions translate between plotly's trace-level axis references and its
layout keys, and read paper coordinates back out of a figure. They are the part
of the plotly backend that has no drawing in it.
"""

import plotly.graph_objects as go
import pytest
from plotly.subplots import make_subplots

from pylocuszoom.backends import plotly_layout


class TestSecondaryAxisKey:
    """Trace-level axis references map to their layout keys."""

    @pytest.mark.parametrize(
        ("secondary_ref", "expected"),
        [("y", "yaxis"), ("y2", "yaxis2"), ("y100", "yaxis100")],
    )
    def test_y_reference_becomes_a_yaxis_key(self, secondary_ref, expected):
        """A y reference gains the 'axis' infix plotly's layout keys carry."""
        assert plotly_layout.secondary_axis_key(secondary_ref) == expected

    def test_a_reference_that_is_not_a_y_axis_passes_through(self):
        """Anything that does not name a y axis is returned unchanged."""
        assert plotly_layout.secondary_axis_key("x2") == "x2"


class TestPanelY:
    """A row's paper coordinates are read from its y-axis domain."""

    @pytest.fixture
    def two_row_figure(self):
        """A two-row subplot figure with explicit, unequal row heights."""
        return make_subplots(rows=2, cols=1, row_heights=[0.75, 0.25])

    @pytest.mark.parametrize("vertical", ["bottom", "middle", "top"])
    def test_coordinate_lies_inside_the_row_domain(self, two_row_figure, vertical):
        """Every anchor of a row falls within that row's own domain."""
        domain = two_row_figure.layout.yaxis.domain

        y = plotly_layout.panel_y(two_row_figure, 1, vertical)

        assert domain[0] <= y <= domain[1]

    def test_the_three_anchors_are_ordered(self, two_row_figure):
        """Bottom, middle and top increase in paper coordinates."""
        anchors = [
            plotly_layout.panel_y(two_row_figure, 1, vertical)
            for vertical in ("bottom", "middle", "top")
        ]

        assert anchors == sorted(anchors)

    def test_middle_is_the_midpoint_of_the_domain(self, two_row_figure):
        """The middle anchor is exactly halfway up the row."""
        domain = two_row_figure.layout.yaxis.domain

        assert plotly_layout.panel_y(two_row_figure, 1, "middle") == pytest.approx(
            (domain[0] + domain[1]) / 2
        )

    def test_a_figure_without_domains_falls_back(self):
        """A plain figure with no subplot domains still yields a coordinate."""
        assert plotly_layout.panel_y(go.Figure(), 1, "bottom") == pytest.approx(0.01)


class TestXRange:
    """A panel's x-range comes from the axis, or from the traces it holds."""

    def test_an_explicit_axis_range_is_used(self):
        """A range set on the axis wins over the data extent."""
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scatter(x=[1, 2, 3], y=[1, 1, 1]), row=1, col=1)
        fig.update_xaxes(range=[0, 10], row=1, col=1)

        panel = plotly_layout._Panel(fig, 1)

        assert tuple(plotly_layout.x_range(panel, "xaxis")) == (0, 10)

    def test_the_trace_extent_is_the_fallback(self):
        """Without an axis range the panel spans its own data."""
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scatter(x=[5, 1, 3], y=[1, 1, 1]), row=1, col=1)

        panel = plotly_layout._Panel(fig, 1)

        assert plotly_layout.x_range(panel, "xaxis") == (1, 5)

    def test_a_panel_with_no_x_data_has_no_range(self):
        """An empty panel reports no range rather than an invented one."""
        panel = plotly_layout._Panel(make_subplots(rows=1, cols=1), 1)

        assert plotly_layout.x_range(panel, "xaxis") is None

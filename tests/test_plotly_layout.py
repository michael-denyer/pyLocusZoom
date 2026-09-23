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


class TestPanelTop:
    """A row's top edge is read from its y-axis domain."""

    def test_each_row_reports_its_own_top(self):
        """The lower row's top is below the upper row's."""
        fig = make_subplots(rows=2, cols=1, row_heights=[0.75, 0.25])

        tops = [
            plotly_layout.panel_top(plotly_layout._Panel(fig, row)) for row in (1, 2)
        ]

        assert tops == [fig.layout.yaxis.domain[1], fig.layout.yaxis2.domain[1]]
        assert tops[1] < tops[0]

    def test_a_figure_without_domains_falls_back(self):
        """A plain figure with no subplot domains still yields a coordinate."""
        assert plotly_layout.panel_top(
            plotly_layout._Panel(go.Figure(), 1)
        ) == pytest.approx(0.99)


class TestXRange:
    """A panel's x-range is its axis limits, shared across a matched group."""

    def test_an_explicit_axis_range_is_used(self):
        """A range set on the axis wins over the data extent."""
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scatter(x=[1, 2, 3], y=[1, 1, 1]), row=1, col=1)
        fig.update_xaxes(range=[0, 10], row=1, col=1)

        panel = plotly_layout._Panel(fig, 1)

        assert tuple(plotly_layout.x_range(panel, "xaxis")) == (0, 10)

    def test_a_matched_axis_takes_the_limits_of_its_group(self):
        """A shared-x panel without limits of its own spans the shared range."""
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
        fig.add_trace(go.Scatter(x=[4, 5], y=[1, 1]), row=1, col=1)
        fig.update_xaxes(range=[0, 10], row=1, col=1)

        lower = plotly_layout._Panel(fig, 2)

        assert tuple(plotly_layout.x_range(lower, "xaxis2")) == (0, 10)

    def test_the_traces_are_never_consulted(self):
        """Without limits anywhere the panel has no range, whatever it holds."""
        fig = make_subplots(rows=1, cols=1)
        fig.add_trace(go.Scatter(x=[5, 1, 3], y=[1, 1, 1]), row=1, col=1)

        panel = plotly_layout._Panel(fig, 1)

        assert plotly_layout.x_range(panel, "xaxis") is None

    def test_a_panel_with_no_x_data_has_no_range(self):
        """An empty panel reports no range rather than an invented one."""
        panel = plotly_layout._Panel(make_subplots(rows=1, cols=1), 1)

        assert plotly_layout.x_range(panel, "xaxis") is None

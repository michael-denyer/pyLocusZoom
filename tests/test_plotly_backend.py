"""Tests for the Plotly backend's subplot addressing and titles."""

import pandas as pd

from pylocuszoom.backends.composition import LegendEntry, mb_tick_positions
from pylocuszoom.backends.plotly_backend import PlotlyBackend
from pylocuszoom.manhattan_plotter import ManhattanPlotter


class TestPlotlyGridSubplotAxisAddressing:
    """Critical: Plotly grid subplots are misaddressed.

    Bug: axis helpers use row-only axis names and several helpers hard-code col=1,
    so in create_figure_grid the QQ column doesn't receive axis limits/labels and
    Manhattan can be overwritten by QQ settings; lines/shapes land only in column 1.
    """

    def test_plotly_axis_name_accounts_for_column(self):
        """_axis_name should return different names for different columns in a grid."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure_grid(n_rows=1, n_cols=2, figsize=(12, 6))

        # In a 1x2 grid:
        # - (row=1, col=1) should use xaxis/yaxis (subplot index 1)
        # - (row=1, col=2) should use xaxis2/yaxis2 (subplot index 2)

        # Current bug: _axis_name only considers row, not column
        # For row=1 it always returns "xaxis"/"yaxis" regardless of column

        manhattan_ax = axes[0]  # (fig, row=1, col=1)
        qq_ax = axes[1]  # (fig, row=1, col=2)

        # Set different y-axis labels
        backend.set_ylabel(manhattan_ax, "Manhattan Y")
        backend.set_ylabel(qq_ax, "QQ Y")

        # Verify they are different axes - check layout has both yaxis and yaxis2
        # If bug exists, both labels go to yaxis and yaxis2 is never set
        layout = fig.layout

        # Check that we have distinct y-axis configurations
        yaxis_title = layout.yaxis.title.text if layout.yaxis.title else None
        yaxis2_title = (
            layout.yaxis2.title.text
            if hasattr(layout, "yaxis2") and layout.yaxis2 and layout.yaxis2.title
            else None
        )

        assert yaxis2_title is not None, (
            "yaxis2 should have a title set for column 2, but it's None. "
            "Bug: _axis_name doesn't account for column."
        )
        assert yaxis_title != yaxis2_title, (
            f"yaxis and yaxis2 have same title '{yaxis_title}'. "
            "Bug: both columns writing to same axis."
        )

    @staticmethod
    def _qq_column_with_a_trace(backend):
        """A 1x2 grid whose column-2 panel holds a point.

        Plotly drops a row/col-addressed hline or vline on a subplot with no
        trace, so the column needs data before a line can land on it.
        """
        fig, axes = backend.create_figure_grid(n_rows=1, n_cols=2, figsize=(12, 6))
        backend.scatter(axes[1], pd.Series([1.0]), pd.Series([1.0]), colors="red")
        return fig, axes[1]

    def test_plotly_axhline_targets_correct_column(self):
        """axhline spans column 2's x domain at a y in column 2's data units."""
        backend = PlotlyBackend()
        fig, qq_ax = self._qq_column_with_a_trace(backend)

        backend.axhline(qq_ax, y=5.0, color="red")

        assert [(s.type, s.xref, s.yref, s.y0) for s in fig.layout.shapes] == [
            ("line", "x2 domain", "y2", 5.0)
        ]

    def test_plotly_add_rectangle_targets_correct_column(self):
        """add_rectangle should add shape to the correct column."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure_grid(n_rows=1, n_cols=2, figsize=(12, 6))

        backend.add_rectangle(axes[1], xy=(0, 0), width=1, height=1)

        assert [(s.type, s.xref, s.yref) for s in fig.layout.shapes] == [
            ("rect", "x2", "y2")
        ]

    def test_plotly_add_polygon_targets_correct_column(self):
        """add_polygon should add shape to the correct column."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure_grid(n_rows=1, n_cols=2, figsize=(12, 6))

        backend.add_polygon(axes[1], points=[[0, 0], [1, 0], [0.5, 1]])

        assert [(s.type, s.xref, s.yref) for s in fig.layout.shapes] == [
            ("path", "x2", "y2")
        ]

    def test_plotly_axvline_targets_correct_column(self):
        """axvline spans column 2's y domain at an x in column 2's data units."""
        backend = PlotlyBackend()
        fig, qq_ax = self._qq_column_with_a_trace(backend)

        backend.axvline(qq_ax, x=5.0, color="red")

        assert [(s.type, s.xref, s.yref, s.x0) for s in fig.layout.shapes] == [
            ("line", "x2", "y2 domain", 5.0)
        ]

    def test_plotly_set_xlim_targets_correct_column(self):
        """set_xlim should set limits on the correct column's x-axis."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure_grid(n_rows=1, n_cols=2, figsize=(12, 6))

        manhattan_ax = axes[0]  # (fig, row=1, col=1)
        qq_ax = axes[1]  # (fig, row=1, col=2)

        # Set different x-limits for each column
        backend.set_xlim(manhattan_ax, 0, 1000)
        backend.set_xlim(qq_ax, 0, 10)

        # Bug: set_xlim uses _axis_name which only considers row
        # Both columns should have different x-axis ranges
        layout = fig.layout

        xaxis_range = layout.xaxis.range if layout.xaxis.range else None
        xaxis2_range = (
            layout.xaxis2.range if hasattr(layout, "xaxis2") and layout.xaxis2 else None
        )

        assert xaxis2_range is not None, (
            "xaxis2 range should be set for column 2, but it's None"
        )
        assert xaxis_range != xaxis2_range, (
            "xaxis and xaxis2 have same range. Bug: both columns using same axis."
        )

    def test_plot_manhattan_qq_distinct_axes(self, manhattan_rs_gwas_df):
        """plot_manhattan_qq should have distinct axis limits for Manhattan and QQ."""
        plotter = ManhattanPlotter(species="canine", backend="plotly")
        fig = plotter.plot_manhattan_qq(manhattan_rs_gwas_df)

        layout = fig.layout

        # Manhattan (column 1) x-axis should have large cumulative positions
        # QQ (column 2) x-axis should have small expected -log10(p) values

        # Check that both axes exist and have different ranges
        xaxis_range = layout.xaxis.range if layout.xaxis.range else None
        xaxis2_range = (
            layout.xaxis2.range
            if hasattr(layout, "xaxis2") and layout.xaxis2 and layout.xaxis2.range
            else None
        )

        assert xaxis_range is not None, "Manhattan x-axis should have range set"
        assert xaxis2_range is not None, "QQ plot x-axis (xaxis2) should have range set"

        # Manhattan range (in bp) should be much larger than QQ range (in -log10(p))
        manhattan_span = xaxis_range[1] - xaxis_range[0]
        qq_span = xaxis2_range[1] - xaxis2_range[0]

        # Manhattan positions are in millions, QQ is typically 0-10
        assert manhattan_span > 1000, (
            f"Manhattan x-range ({manhattan_span}) should be large (genomic positions)"
        )
        assert qq_span < 100, (
            f"QQ x-range ({qq_span}) should be small (-log10(p) values)"
        )

    @staticmethod
    def _grid_2x2():
        backend = PlotlyBackend()
        fig, axes = backend.create_figure_grid(n_rows=2, n_cols=2)
        for ax in axes:  # add_vrect skips subplots that hold no trace
            backend.scatter(ax, pd.Series([0.0, 30.0]), pd.Series([0.0, 1.0]), "blue")
        return backend, fig, axes

    def test_region_highlight_lands_on_the_panels_given(self):
        backend, fig, axes = self._grid_2x2()

        backend.add_region_highlight([axes[1], axes[3]], 10, 20)

        assert [(s.xref, s.yref) for s in fig.layout.shapes] == [
            ("x2", "y2 domain"),
            ("x4", "y4 domain"),
        ]

    def test_legend_lands_on_the_panel_given(self):
        backend, fig, axes = self._grid_2x2()

        backend.add_legend(axes[3], [LegendEntry(label="a", color="red", marker="o")])

        entry = fig.data[-1]
        assert (entry.xaxis, entry.yaxis) == ("x4", "y4")
        assert fig.layout.legend.y == fig.layout.yaxis4.domain[1]

    def test_fill_between_accepts_a_scalar_y2(self):
        backend = PlotlyBackend()
        fig, axes = backend.create_figure(height_ratios=[1.0], figsize=(6, 4))

        backend.fill_between(axes[0], pd.Series([1.0, 2.0]), pd.Series([3.0, 4.0]), 0.5)

        assert list(fig.data[-1].y) == [0.5, 0.5, 4.0, 3.0]


class TestPlotlySetTitleOverwriting:
    """Low: Plotly set_title only updates the overall figure title for row 1.

    Bug: In plot_manhattan_qq the Manhattan title is overwritten by the QQ title,
    and in stacked mode only the first row gets a QQ title.
    """

    def test_plotly_set_title_per_subplot(self):
        """set_title should set title for specific subplot using annotations for grids."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure_grid(n_rows=1, n_cols=2, figsize=(12, 6))

        manhattan_ax = axes[0]
        qq_ax = axes[1]

        # Set titles for each subplot
        backend.set_title(manhattan_ax, "Manhattan Plot")
        backend.set_title(qq_ax, "QQ Plot")

        # For grid layouts, titles are now added as annotations
        annotations = fig.layout.annotations
        assert annotations is not None and len(annotations) >= 2, (
            f"Expected at least 2 annotations for subplot titles, got {len(annotations) if annotations else 0}"
        )

        # Extract annotation texts
        annotation_texts = [ann.text for ann in annotations]

        # Both titles should appear (with potential HTML formatting)
        assert any("Manhattan" in str(t) for t in annotation_texts), (
            f"Manhattan title not found in annotations: {annotation_texts}"
        )
        assert any("QQ" in str(t) for t in annotation_texts), (
            f"QQ title not found in annotations: {annotation_texts}"
        )

    def test_plot_manhattan_qq_has_distinct_titles(self, manhattan_rs_gwas_df):
        """plot_manhattan_qq should show both Manhattan and QQ titles."""
        plotter = ManhattanPlotter(species="canine", backend="plotly")
        fig = plotter.plot_manhattan_qq(manhattan_rs_gwas_df)

        # Convert to JSON to inspect all text elements
        import json

        fig_json = json.loads(fig.to_json())

        # Look for title text in annotations
        all_text = []

        # Check annotations (grid layouts use annotations for titles)
        for ann in fig_json.get("layout", {}).get("annotations", []):
            if "text" in ann:
                all_text.append(ann["text"])

        # We should see both plot types in annotations
        text_combined = " ".join(all_text).lower()

        has_manhattan = "manhattan" in text_combined
        has_qq = "qq" in text_combined or "λ" in text_combined

        assert has_manhattan, f"Manhattan title not found in annotations: {all_text}"
        assert has_qq, f"QQ title not found in annotations: {all_text}"


class TestPlotlyMegabaseTicksFollowTheAxisRange:
    """Tick ladders come from the axis limits, never from a panel's traces."""

    def test_each_panel_is_ticked_from_its_own_limits(self):
        backend = PlotlyBackend()
        fig, panels = backend.create_figure(
            height_ratios=[1.0, 1.0], figsize=(8, 6), sharex=False
        )
        backend.set_xlim(panels[0], 1_000_000, 2_000_000)
        backend.set_xlim(panels[1], 50_000_000, 51_000_000)

        backend.format_xaxis_mb(panels[0])
        backend.format_xaxis_mb(panels[1])

        assert (
            list(fig.layout.xaxis.tickvals)
            == mb_tick_positions(1_000_000, 2_000_000)[0]
        )
        assert (
            list(fig.layout.xaxis2.tickvals)
            == mb_tick_positions(50_000_000, 51_000_000)[0]
        )

    def test_shared_axis_panel_without_its_own_limits_is_ticked_from_the_shared_range(
        self,
    ):
        """A stacked panel whose traces are narrower than the region must still
        carry the region's tick ladder, taken from the axis it is matched to."""
        backend = PlotlyBackend()
        fig, panels = backend.create_figure(
            height_ratios=[1.0, 1.0], figsize=(8, 6), sharex=True
        )
        backend.set_xlim(panels[0], 1_000_000, 2_000_000)
        backend.scatter(
            panels[1], pd.Series([1_410_000, 1_660_000]), pd.Series([1.0, 2.0]), "blue"
        )

        backend.format_xaxis_mb(panels[0])
        backend.format_xaxis_mb(panels[1])

        expected = mb_tick_positions(1_000_000, 2_000_000)[0]
        assert list(fig.layout.xaxis.tickvals) == expected
        assert list(fig.layout.xaxis2.tickvals) == expected

    def test_panel_with_no_range_anywhere_is_left_to_plotly(self):
        backend = PlotlyBackend()
        fig, panels = backend.create_figure(height_ratios=[1.0], figsize=(8, 6))
        backend.scatter(
            panels[0], pd.Series([1_410_000, 1_660_000]), pd.Series([1.0, 2.0]), "blue"
        )

        backend.format_xaxis_mb(panels[0])

        assert fig.layout.xaxis.tickvals is None

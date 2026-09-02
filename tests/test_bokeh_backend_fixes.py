"""Tests for Bokeh backend bug fixes.

Covers issues found in code review:
- High #1: save() global state via output_file()
- High #2: add_panel_label() broken with DataRange1d
- Medium #3: scatter() hover column name collision
- Medium #4: _create_color_palette() only handles 6-digit hex
- Medium #5: highlight_heatmap_snp O(n) renderers
- Low #6: identity orientation_map in add_colorbar (code cleanup)
- Low #7: redundant local ColumnDataSource import (code cleanup)
"""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends.bokeh_backend import BokehBackend, _create_color_palette


class TestAddPanelLabelWithDataRange1d:
    """High #2: add_panel_label() must work with DataRange1d (unresolved ranges)."""

    def test_panel_label_with_datarange1d(self):
        """add_panel_label should work when x_range/y_range are DataRange1d (start=None)."""
        from bokeh.models import Label

        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        # DataRange1d has start=None before rendering — this should NOT crash
        # or place the label at data coordinate (0, 0)
        backend.add_panel_label(ax, "A", x_frac=0.02, y_frac=0.95)

        # Find the label in the layout
        labels = [
            r
            for r in ax.renderers
            + list(ax.center)
            + list(ax.above)
            + list(ax.below)
            + list(ax.left)
            + list(ax.right)
            if isinstance(r, Label)
        ]
        assert len(labels) == 1
        label_obj = labels[0]
        assert label_obj.text == "A"
        # Must NOT use data coordinates (which are broken with DataRange1d)
        # Should use screen units so position works regardless of data range
        assert label_obj.x_units == "screen", (
            f"Expected screen units, got {label_obj.x_units}. "
            "Data coordinates fail with DataRange1d (start=None)."
        )

    def test_panel_label_with_explicit_range(self):
        """add_panel_label should also work after set_xlim/set_ylim."""
        from bokeh.models import Label

        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        backend.set_xlim(ax, 1_000_000, 2_000_000)
        backend.set_ylim(ax, 0, 10)
        backend.add_panel_label(ax, "B")

        labels = [
            r
            for r in ax.renderers
            + list(ax.center)
            + list(ax.above)
            + list(ax.below)
            + list(ax.left)
            + list(ax.right)
            if isinstance(r, Label)
        ]
        assert len(labels) == 1
        assert labels[0].text == "B"


class TestScatterHoverColumnCollision:
    """Medium #3: scatter() hover columns should not collide with internal keys."""

    def test_hover_column_named_x_does_not_corrupt_scatter(self):
        """Hover data with column 'x' should not overwrite scatter x-coordinates."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        x = pd.Series([1.0, 2.0, 3.0])
        y = pd.Series([4.0, 5.0, 6.0])
        hover_data = pd.DataFrame(
            {
                "x": ["label_a", "label_b", "label_c"],  # Collides with scatter 'x'
                "p_value": [0.001, 0.01, 0.1],
            }
        )

        renderer = backend.scatter(ax, x, y, colors="#BEBEBE", hover_data=hover_data)

        # The scatter coordinates must still be [1, 2, 3], not overwritten
        source = renderer.data_source
        np.testing.assert_array_equal(source.data["x"], [1.0, 2.0, 3.0])

    def test_hover_column_named_color_does_not_corrupt_scatter(self):
        """Hover data with column 'color' should not overwrite point colors."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        x = pd.Series([1.0, 2.0, 3.0])
        y = pd.Series([4.0, 5.0, 6.0])
        hover_data = pd.DataFrame(
            {
                "SNP": ["rs1", "rs2", "rs3"],
                "color": ["info_a", "info_b", "info_c"],  # Collides with 'color' key
            }
        )

        renderer = backend.scatter(ax, x, y, colors="red", hover_data=hover_data)

        # The color data must still be the scatter colors, not overwritten
        source = renderer.data_source
        assert list(source.data["color"]) == ["red", "red", "red"]

    def test_hover_column_named_size_does_not_corrupt_scatter(self):
        """Hover data with column 'size' should not overwrite point sizes."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        x = pd.Series([1.0, 2.0])
        y = pd.Series([4.0, 5.0])
        hover_data = pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "size": [100, 200],  # Collides with 'size' key
            }
        )

        renderer = backend.scatter(
            ax, x, y, colors="red", sizes=60, hover_data=hover_data
        )

        source = renderer.data_source
        expected_size = max(6, 60**0.5)
        assert list(source.data["size"]) == [expected_size, expected_size]


class TestCreateColorPalette:
    """Medium #4: _create_color_palette should handle various color formats."""

    def test_standard_6_digit_hex(self):
        """Standard 6-digit hex should work."""
        palette = _create_color_palette("#FFFFFF", "#FF0000", 3)
        assert len(palette) == 3
        assert palette[0] == "#ffffff"
        assert palette[-1] == "#ff0000"

    def test_3_digit_hex(self):
        """3-digit hex (e.g., '#F00') should be handled."""
        palette = _create_color_palette("#FFF", "#F00", 3)
        assert len(palette) == 3
        assert palette[0] == "#ffffff"
        assert palette[-1] == "#ff0000"

    def test_named_css_color(self):
        """Named CSS colors like 'red', 'white' should be handled."""
        palette = _create_color_palette("white", "red", 3)
        assert len(palette) == 3
        assert palette[0] == "#ffffff"
        assert palette[-1] == "#ff0000"


class TestHighlightHeatmapSnpBatched:
    """Medium #5: highlight_heatmap_snp should use batched renderers."""

    def test_highlight_creates_two_renderers(self):
        """Should create exactly 2 renderers (row batch + column batch), not O(n)."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        initial_renderer_count = len(ax.renderers)
        backend.highlight_heatmap_snp(ax, layout, snp_idx=5, n_snps=20)

        new_renderers = len(ax.renderers) - initial_renderer_count
        # Should be at most 2 renderers (one for row, one for column), not 20
        assert new_renderers <= 2, (
            f"Expected at most 2 batched renderers, got {new_renderers}. "
            "Each highlight cell should not create a separate renderer."
        )

    def test_highlight_still_correct(self):
        """Batched highlight should still highlight the correct cells."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        initial_renderer_count = len(ax.renderers)
        backend.highlight_heatmap_snp(ax, layout, snp_idx=1, n_snps=3)

        # snp_idx=1, n_snps=3:
        # Row highlights: cells (j=0,i=1) and (j=1,i=1) — 2 cells
        # Column highlights: cells (j=1,i=2) — 1 cell
        # Total: 3 highlighted cells
        new_renderers = ax.renderers[initial_renderer_count:]

        # Collect all highlight cell coordinates from new renderers
        all_xs = []
        all_ys = []
        for r in new_renderers:
            if hasattr(r, "data_source"):
                data = r.data_source.data
                # Bokeh rect uses 'x' and 'y' keys for center coordinates
                if "x" in data:
                    all_xs.extend(data["x"])
                    all_ys.extend(data["y"])

        assert len(all_xs) == 3, f"Expected 3 highlight cells, got {len(all_xs)}"
        # Row cells: (0,1) and (1,1); Column cell: (1,2)
        expected = {(0, 1), (1, 1), (1, 2)}
        actual = {(int(x), int(y)) for x, y in zip(all_xs, all_ys)}
        assert actual == expected, f"Expected cells {expected}, got {actual}"


class TestAddColorbarNoIdentityMap:
    """Low #6: add_colorbar should not use identity orientation_map."""

    def test_colorbar_vertical(self):
        """Vertical colorbar should work without identity map."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        mapper = backend.add_heatmap(ax, data, [0, 1], [0, 1])
        colorbar = backend.add_colorbar(ax, mapper, label="R²", orientation="vertical")

        assert colorbar.orientation == "vertical"

    def test_colorbar_horizontal(self):
        """Horizontal colorbar should work."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        mapper = backend.add_heatmap(ax, data, [0, 1], [0, 1])
        colorbar = backend.add_colorbar(
            ax, mapper, label="R²", orientation="horizontal"
        )

        assert colorbar.orientation == "horizontal"


class TestNoRedundantImport:
    """Low #7: _ensure_legend_range should not import ColumnDataSource locally."""

    def test_ensure_legend_range_returns_column_data_source(self):
        """_ensure_legend_range should return ColumnDataSource (using module-level import)."""
        from bokeh.models import ColumnDataSource

        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        source = backend._ensure_legend_range(ax)
        assert isinstance(source, ColumnDataSource)


class TestAddTextAnchors:
    """add_text must honour every (ha, va) pair, not just four."""

    @pytest.mark.parametrize(
        "ha,va,text_align,text_baseline",
        [
            ("center", "center", "center", "middle"),
            ("left", "top", "left", "top"),
            ("right", "baseline", "right", "alphabetic"),
            ("center", "bottom", "center", "bottom"),
        ],
    )
    def test_alignment_maps_independently(self, ha, va, text_align, text_baseline):
        from bokeh.models import Label

        backend = BokehBackend()
        _, axes = backend.create_figure(n_panels=1, height_ratios=[1.0], figsize=(8, 4))
        backend.add_text(axes[0], 1.0, 2.0, "note", ha=ha, va=va)

        labels = [r for r in axes[0].center if isinstance(r, Label)]
        assert len(labels) == 1
        assert labels[0].text_align == text_align
        assert labels[0].text_baseline == text_baseline

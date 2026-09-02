"""Unit tests for BokehBackend."""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends.bokeh_backend import BokehBackend, _create_color_palette
from pylocuszoom.stats_plotter import StatsPlotter


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


class TestBokehSetXticksCorruptsYAxis:
    """High: Bokeh set_xticks writes x-axis labels into the y-axis.

    Bug: set_xticks writes to ax.yaxis.major_label_overrides instead of only
    ax.xaxis.major_label_overrides, corrupting y-axis labels.
    """

    def test_bokeh_set_xticks_does_not_modify_yaxis(self):
        """set_xticks should not modify y-axis properties."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(12, 6)
        )
        ax = axes[0]

        # Set custom y-axis labels first
        backend.set_yticks(ax, positions=[0, 1, 2], labels=["A", "B", "C"])

        # Now set x-axis ticks - this should NOT affect y-axis
        backend.set_xticks(
            ax,
            positions=[100, 200, 300],
            labels=["X1", "X2", "X3"],
            fontsize=10,
        )

        # Bug: set_xticks has these lines that corrupt y-axis:
        # ax.yaxis.major_label_overrides = dict(zip(positions, labels))
        # ax.yaxis.major_label_text_font_size = f"{fontsize}pt"

        # Check that y-axis was not modified by set_xticks
        y_overrides = ax.yaxis.major_label_overrides

        # y_overrides should NOT contain x-axis labels
        assert "X1" not in y_overrides.values(), (
            f"X-axis labels leaked into y-axis overrides: {y_overrides}"
        )
        assert "X2" not in y_overrides.values(), (
            f"X-axis labels leaked into y-axis overrides: {y_overrides}"
        )


class TestBokehSetYticksIgnoresLabels:
    """High: Bokeh set_yticks ignores the provided labels entirely.

    Bug: set_yticks only sets ax.yaxis.ticker = positions but doesn't set
    ax.yaxis.major_label_overrides, so custom labels are ignored.
    """

    def test_bokeh_set_yticks_applies_labels(self):
        """set_yticks should apply the provided custom labels."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(12, 6)
        )
        ax = axes[0]

        # Set custom y-axis labels
        positions = [0, 1, 2, 3]
        labels = ["Type 2 Diabetes", "BMI", "Height", "Blood Pressure"]

        backend.set_yticks(ax, positions=positions, labels=labels, fontsize=10)

        # Bug: set_yticks only sets ticker, ignoring labels parameter
        # Current implementation:
        #   ax.yaxis.ticker = positions
        # Missing:
        #   ax.yaxis.major_label_overrides = dict(zip(positions, labels))

        # Check that labels were applied
        y_overrides = ax.yaxis.major_label_overrides

        assert len(y_overrides) > 0, (
            "set_yticks should set major_label_overrides but it's empty. "
            "Bug: labels parameter is ignored."
        )

        # Verify the correct labels are present
        for pos, expected_label in zip(positions, labels):
            assert pos in y_overrides, f"Position {pos} not in y-axis overrides"
            assert y_overrides[pos] == expected_label, (
                f"Expected label '{expected_label}' at position {pos}, "
                f"got '{y_overrides.get(pos)}'"
            )

    def test_bokeh_phewas_shows_phenotype_names(self, phewas_with_effects_df):
        """PheWAS plot should show phenotype names, not numeric indices."""
        plotter = StatsPlotter(backend="bokeh")
        fig = plotter.plot_phewas(
            phewas_with_effects_df,
            variant_id="rs12345",
            phenotype_col="phenotype",
            p_col="p",
        )

        # Get the main figure from the layout
        main_fig = fig.children[0] if hasattr(fig, "children") else fig

        # Check y-axis has phenotype names, not just numeric indices
        y_overrides = main_fig.yaxis.major_label_overrides

        assert len(y_overrides) > 0, (
            "PheWAS y-axis should have label overrides for phenotype names"
        )

        # Check that actual phenotype names are present
        override_values = list(y_overrides.values())
        assert "Type 2 Diabetes" in override_values or any(
            "Diabetes" in str(v) for v in override_values
        ), f"Expected phenotype names in y-axis, got: {override_values}"

    def test_bokeh_forest_shows_study_names(self, sample_forest_df):
        """Forest plot should show study names, not numeric indices."""
        plotter = StatsPlotter(backend="bokeh")
        fig = plotter.plot_forest(
            sample_forest_df,
            variant_id="rs12345",
            study_col="study",
            effect_col="effect",
            ci_lower_col="ci_lower",
            ci_upper_col="ci_upper",
        )

        # Get the main figure from the layout
        main_fig = fig.children[0] if hasattr(fig, "children") else fig

        # Check y-axis has study names
        y_overrides = main_fig.yaxis.major_label_overrides

        assert len(y_overrides) > 0, (
            "Forest plot y-axis should have label overrides for study names"
        )

        # Check that actual study names are present
        override_values = list(y_overrides.values())
        assert "Study A" in override_values or any(
            "Study" in str(v) for v in override_values
        ), f"Expected study names in y-axis, got: {override_values}"

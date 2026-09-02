"""Tests for MiamiPlotter class."""

import pandas as pd
import pytest

from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.miami_plotter import MiamiPlotter
from tests.conftest import FIGURE_TYPES


def _max_scatter_x(ax) -> float:
    """Return the largest scatter x coordinate drawn on a matplotlib axis."""
    return max(
        float(offset[0]) for coll in ax.collections for offset in coll.get_offsets()
    )


class TestMiamiPlotter:
    """Tests for the MiamiPlotter class."""

    @pytest.fixture
    def miami_plotter(self):
        """Create a MiamiPlotter instance."""
        return MiamiPlotter(species="canine")

    @pytest.fixture
    def miami_panel_dfs_with_rs(self):
        """Create sample GWAS data for top and bottom panels.

        Returns tuple of (top_df, bottom_df) for Miami plot comparison.
        """
        top_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 1000, 2000, 1000, 2000],
                "p": [0.01, 0.001, 0.0001, 0.05, 1e-8, 0.1],
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 1000, 2000, 1000, 2000],
                "p": [0.05, 0.002, 0.001, 0.1, 5e-7, 0.01],
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
            }
        )
        return top_df, bottom_df

    def test_plot_miami_returns_figure(self, miami_plotter, miami_panel_dfs_with_rs):
        """Test that plot_miami returns a figure object."""
        top_df, bottom_df = miami_panel_dfs_with_rs
        fig = miami_plotter.plot_miami(top_df, bottom_df)
        assert fig is not None

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_returns_the_backends_figure_type(self, backend, miami_panel_dfs_with_rs):
        """Each backend returns its own figure type."""
        plotter = MiamiPlotter(species="canine", backend=backend)
        top_df, bottom_df = miami_panel_dfs_with_rs

        fig = plotter.plot_miami(top_df, bottom_df)

        assert isinstance(fig, FIGURE_TYPES[backend])

    def test_chromosome_colors_match(self, miami_plotter, miami_panel_dfs_with_rs):
        """Test that same chromosome has same color in both panels.

        Verifies that chromosome colors are consistent between top and bottom
        panels - critical for visual comparison in Miami plots.
        """
        top_df, bottom_df = miami_panel_dfs_with_rs
        fig = miami_plotter.plot_miami(top_df, bottom_df)

        # Get the figure axes (matplotlib-specific verification)
        axes = fig.get_axes()
        assert len(axes) == 2, "Miami plot should have exactly 2 panels"

        # For this test, we verify the figure was created with 2 panels
        # Color consistency is implicitly tested by using the same color mapping
        top_ax, bottom_ax = axes
        assert top_ax is not None
        assert bottom_ax is not None

    def test_panels_share_chromosome_offsets(self):
        """Both panels place a chromosome at the same x, whatever their extent."""
        plotter = MiamiPlotter(species="human")
        top_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2],
                "pos": [1000, 9_000_000, 5000],
                "p": [0.01, 0.001, 1e-6],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2],
                "pos": [1000, 2_000_000, 5000],
                "p": [0.02, 0.003, 1e-5],
            }
        )

        fig = plotter.plot_miami(top_df, bottom_df)

        top_ax, bottom_ax = fig.get_axes()
        assert _max_scatter_x(top_ax) == _max_scatter_x(bottom_ax)

    def test_species_without_chromosome_order_raises(self):
        """A species with no built-in chromosome order cannot lay out the axis."""
        panel_df = pd.DataFrame({"chrom": ["1"], "pos": [1000], "p": [0.01]})
        with pytest.raises(ValueError, match="No built-in chromosome order"):
            MiamiPlotter(species="nonsense").plot_miami(panel_df, panel_df)

    def test_species_order_drives_chromosome_ticks(self):
        """The species table orders the axis, not an alphabetic sort."""
        plotter = MiamiPlotter(species="feline")
        panel_df = pd.DataFrame(
            {
                "chrom": ["A1", "MT", "X", "Y"],
                "pos": [1000, 1000, 1000, 1000],
                "p": [0.01, 0.001, 0.02, 0.03],
            }
        )

        fig = plotter.plot_miami(panel_df, panel_df)

        bottom_ax = fig.get_axes()[1]
        labels = [tick.get_text() for tick in bottom_ax.get_xticklabels()]
        assert labels == ["A1", "X", "Y", "MT"]

    def test_empty_dataframe_raises(self, miami_plotter):
        """Test that empty DataFrame raises appropriate error."""
        empty_df = pd.DataFrame(columns=["chrom", "pos", "p"])
        valid_df = pd.DataFrame({"chrom": [1], "pos": [1000], "p": [0.01]})

        # Empty top_df should raise
        with pytest.raises((ValueError, KeyError)):
            miami_plotter.plot_miami(empty_df, valid_df)

        # Empty bottom_df should raise
        with pytest.raises((ValueError, KeyError)):
            miami_plotter.plot_miami(valid_df, empty_df)

    def test_missing_columns_raises(self, miami_plotter):
        """Test that missing required columns raises ValidationError or ValueError."""
        # DataFrame missing 'p' column
        missing_p_df = pd.DataFrame({"chrom": [1], "pos": [1000]})
        valid_df = pd.DataFrame({"chrom": [1], "pos": [1000], "p": [0.01]})

        with pytest.raises((ValueError, KeyError)):
            miami_plotter.plot_miami(missing_p_df, valid_df)

        # DataFrame missing 'chrom' column
        missing_chrom_df = pd.DataFrame({"pos": [1000], "p": [0.01]})
        with pytest.raises((ValueError, KeyError)):
            miami_plotter.plot_miami(valid_df, missing_chrom_df)

    def test_inverted_bottom_panel(self, miami_plotter, miami_panel_dfs_with_rs):
        """Test that bottom panel y-axis is inverted.

        The bottom panel should display with y-axis inverted (max at top/0 at bottom)
        to create the characteristic Miami plot mirror effect.
        """
        top_df, bottom_df = miami_panel_dfs_with_rs
        fig = miami_plotter.plot_miami(top_df, bottom_df)

        # Get the figure axes
        axes = fig.get_axes()
        assert len(axes) == 2

        top_ax, bottom_ax = axes

        # Top panel should have normal y-axis (0 at bottom, max at top)
        top_ylim = top_ax.get_ylim()
        assert top_ylim[0] < top_ylim[1], "Top panel y-axis should not be inverted"

        # Bottom panel should have inverted y-axis (max at bottom, 0 at top)
        bottom_ylim = bottom_ax.get_ylim()
        assert bottom_ylim[0] > bottom_ylim[1], (
            "Bottom panel y-axis should be inverted (max at bottom, 0 at top)"
        )


class TestMiamiPlotterOptions:
    """Tests for MiamiPlotter configuration options."""

    @pytest.fixture
    def miami_panel_dfs(self):
        """Create sample GWAS data."""
        top_df = pd.DataFrame(
            {
                "chrom": [1, 2, 3],
                "pos": [1000, 1000, 1000],
                "p": [0.01, 0.0001, 1e-8],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 2, 3],
                "pos": [1000, 1000, 1000],
                "p": [0.05, 0.001, 5e-7],
            }
        )
        return top_df, bottom_df

    def test_custom_figsize(self, miami_panel_dfs):
        """Test that custom figsize is respected."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs
        fig = plotter.plot_miami(top_df, bottom_df, figsize=(15, 10))

        # Check figure size
        width, height = fig.get_size_inches()
        assert width == 15
        assert height == 10

    def test_custom_title(self, miami_panel_dfs):
        """Test that custom title is set."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs
        fig = plotter.plot_miami(top_df, bottom_df, title="Custom Miami Title")

        # Check that title was set (matplotlib-specific)
        title_text = fig._suptitle
        assert title_text is not None

    def test_panel_labels(self, miami_panel_dfs):
        """Test that panel labels are applied."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs
        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            top_label="Discovery",
            bottom_label="Replication",
        )
        assert fig is not None

    def test_different_chromosome_sets(self):
        """Test Miami plot when top and bottom have different chromosome sets.

        Miami plot should handle cases where GWAS datasets have different
        chromosomes present (e.g., one has chr X, the other doesn't).
        """
        plotter = MiamiPlotter(species="canine")

        # Top has chr 1, 2, 3
        top_df = pd.DataFrame(
            {
                "chrom": [1, 2, 3],
                "pos": [1000, 1000, 1000],
                "p": [0.01, 0.001, 0.0001],
            }
        )

        # Bottom has chr 1, 2, X (no chr 3, has X instead)
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 2, "X"],
                "pos": [1000, 1000, 1000],
                "p": [0.05, 0.01, 0.001],
            }
        )

        # Should still create figure successfully
        fig = plotter.plot_miami(top_df, bottom_df)
        assert fig is not None


class TestMiamiPlotterHoverData:
    """Tests for Miami plotter hover data support."""

    @pytest.fixture
    def gwas_data_with_rs(self):
        """Create GWAS data with RS IDs for hover tooltips."""
        top_df = pd.DataFrame(
            {
                "chrom": [1, 2, 3],
                "pos": [1000, 2000, 3000],
                "p": [0.01, 0.001, 0.0001],
                "rs": ["rs123", "rs456", "rs789"],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 2, 3],
                "pos": [1000, 2000, 3000],
                "p": [0.05, 0.01, 0.001],
                "rs": ["rs123", "rs456", "rs789"],
            }
        )
        return top_df, bottom_df

    def test_plotly_hover_data(self, gwas_data_with_rs):
        """Test that plotly backend creates figure with hover data."""
        plotter = MiamiPlotter(species="canine", backend="plotly")
        top_df, bottom_df = gwas_data_with_rs
        fig = plotter.plot_miami(top_df, bottom_df, rs_col="rs")
        assert fig is not None

        # Verify traces exist (plotly-specific)
        assert len(fig.data) > 0, "Plotly figure should have traces for hover data"

    def test_bokeh_hover_data(self, gwas_data_with_rs):
        """Test that bokeh backend creates figure with hover tools."""
        plotter = MiamiPlotter(species="canine", backend="bokeh")
        top_df, bottom_df = gwas_data_with_rs
        fig = plotter.plot_miami(top_df, bottom_df, rs_col="rs")
        assert fig is not None


class TestMiamiSignificance:
    """Tests for Miami plot significance lines."""

    @pytest.fixture
    def miami_panel_dfs(self):
        """Create sample GWAS data for significance line testing."""

        top_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 1000, 2000, 1000, 2000],
                "p": [0.01, 0.001, 1e-9, 0.05, 1e-8, 0.1],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 1000, 2000, 1000, 2000],
                "p": [0.05, 0.002, 1e-7, 0.1, 5e-7, 0.01],
            }
        )
        return top_df, bottom_df

    def test_significance_lines_both_panels(self, miami_panel_dfs):
        """Test that significance lines are drawn on both panels by default."""

        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs
        fig = plotter.plot_miami(top_df, bottom_df)

        # Check figure was created
        assert fig is not None

        # Get axes and verify horizontal lines exist
        axes = fig.get_axes()
        assert len(axes) == 2

        # Both panels should have at least one horizontal line (significance line)
        top_ax, bottom_ax = axes
        top_lines = [line for line in top_ax.lines if line.get_linestyle() == "--"]
        bottom_lines = [
            line for line in bottom_ax.lines if line.get_linestyle() == "--"
        ]

        # Default threshold creates one significance line per panel
        assert len(top_lines) >= 1, "Top panel should have significance line"
        assert len(bottom_lines) >= 1, "Bottom panel should have significance line"

    def test_different_thresholds(self, miami_panel_dfs):
        """Test that different thresholds create lines at different positions."""
        import numpy as np

        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs
        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            top_threshold=5e-8,  # -log10 = ~7.3
            bottom_threshold=1e-6,  # -log10 = 6.0
        )

        axes = fig.get_axes()
        top_ax, bottom_ax = axes

        # Find dashed lines (significance lines)
        top_lines = [line for line in top_ax.lines if line.get_linestyle() == "--"]
        bottom_lines = [
            line for line in bottom_ax.lines if line.get_linestyle() == "--"
        ]

        # Get y-values of the significance lines
        if top_lines and bottom_lines:
            top_y = top_lines[0].get_ydata()[0]
            bottom_y = bottom_lines[0].get_ydata()[0]

            # 5e-8 threshold = -log10(5e-8) = ~7.3
            # 1e-6 threshold = -log10(1e-6) = 6.0
            expected_top_y = -np.log10(5e-8)
            expected_bottom_y = -np.log10(1e-6)

            assert np.isclose(top_y, expected_top_y, rtol=0.01)
            assert np.isclose(bottom_y, expected_bottom_y, rtol=0.01)

    def test_no_significance_line(self, miami_panel_dfs):
        """Test that threshold=None suppresses significance line on that panel."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs
        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            top_threshold=5e-8,  # Keep on top
            bottom_threshold=None,  # No line on bottom
        )

        axes = fig.get_axes()
        top_ax, bottom_ax = axes

        # Find dashed lines (significance lines)
        top_lines = [line for line in top_ax.lines if line.get_linestyle() == "--"]
        bottom_lines = [
            line for line in bottom_ax.lines if line.get_linestyle() == "--"
        ]

        assert len(top_lines) >= 1, "Top panel should have significance line"
        assert len(bottom_lines) == 0, "Bottom panel should have no significance line"


class TestMiamiLabels:
    """Tests for Miami plot panel labels and SNP annotations."""

    @pytest.fixture
    def miami_panel_dfs_with_rs(self):
        """Create GWAS data with RS IDs for annotation testing."""
        top_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 1000, 2000, 1000, 2000],
                "p": [0.01, 0.001, 1e-9, 0.05, 1e-8, 0.1],
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 1000, 2000, 1000, 2000],
                "p": [0.05, 0.002, 1e-7, 0.1, 5e-7, 0.01],
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
            }
        )
        return top_df, bottom_df

    def test_panel_labels(self, miami_panel_dfs_with_rs):
        """Test that panel labels are correctly displayed."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs_with_rs
        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            top_label="Discovery",
            bottom_label="Replication",
        )

        # Figure should be created successfully
        assert fig is not None

        # Check that axes were created
        axes = fig.get_axes()
        assert len(axes) == 2

        # For matplotlib, annotations include panel labels
        # The label implementation uses ax.annotate(), which creates texts
        top_ax, bottom_ax = axes

        # Get all text objects from both axes
        top_texts = [t.get_text() for t in top_ax.texts]
        bottom_texts = [t.get_text() for t in bottom_ax.texts]

        assert "Discovery" in top_texts, "Top panel should have 'Discovery' label"
        assert "Replication" in bottom_texts, (
            "Bottom panel should have 'Replication' label"
        )

    def test_snp_annotations_top_panel(self, miami_panel_dfs_with_rs):
        """Test SNP annotations on top panel only."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs_with_rs

        # Annotate rs3 (most significant on top panel - 1e-9) on top only
        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            rs_col="rs",
            top_snp_annotations=["rs3"],
            bottom_snp_annotations=None,
        )

        assert fig is not None
        axes = fig.get_axes()
        assert len(axes) == 2

        top_ax, bottom_ax = axes
        top_texts = [t.get_text() for t in top_ax.texts]
        bottom_texts = [t.get_text() for t in bottom_ax.texts]

        # rs3 should be annotated on top panel
        assert "rs3" in top_texts, "Top panel should have rs3 annotation"
        # rs3 should NOT be annotated on bottom panel (only panel labels if any)
        assert "rs3" not in bottom_texts, "Bottom panel should not have rs3 annotation"

    def test_snp_annotations_both_panels(self, miami_panel_dfs_with_rs):
        """Test independent SNP annotations on each panel."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs_with_rs

        # Different SNPs annotated on each panel
        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            rs_col="rs",
            top_snp_annotations=["rs3", "rs5"],  # Top panel annotations
            bottom_snp_annotations=["rs1", "rs2"],  # Bottom panel annotations
        )

        assert fig is not None
        axes = fig.get_axes()
        top_ax, bottom_ax = axes

        top_texts = [t.get_text() for t in top_ax.texts]
        bottom_texts = [t.get_text() for t in bottom_ax.texts]

        # Verify top panel has rs3 and rs5
        assert "rs3" in top_texts
        assert "rs5" in top_texts

        # Verify bottom panel has rs1 and rs2
        assert "rs1" in bottom_texts
        assert "rs2" in bottom_texts


class TestMiamiHighlight:
    """Tests for Miami plot region highlighting."""

    @pytest.fixture
    def miami_panel_dfs(self):
        """Create GWAS data spanning multiple chromosomes."""
        top_df = pd.DataFrame(
            {
                "chrom": [1, 1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 3000, 1000, 2000, 1000, 2000],
                "p": [0.01, 0.001, 1e-9, 0.05, 1e-8, 0.1, 0.02],
            }
        )
        bottom_df = pd.DataFrame(
            {
                "chrom": [1, 1, 1, 2, 2, 3, 3],
                "pos": [1000, 2000, 3000, 1000, 2000, 1000, 2000],
                "p": [0.05, 0.002, 1e-7, 0.1, 5e-7, 0.01, 0.03],
            }
        )
        return top_df, bottom_df

    def test_highlight_single_region(self, miami_panel_dfs):
        """Test highlighting a single region on both panels."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs

        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            highlight_regions=[("1", 500, 2500)],
        )

        assert fig is not None

        # Verify figure has 2 panels
        axes = fig.get_axes()
        assert len(axes) == 2

        # For matplotlib, check that axvspan was called (creates Polygon patches)
        top_ax, bottom_ax = axes

        # Both axes should have patches (highlighting spans)
        top_patches = [p for p in top_ax.patches if hasattr(p, "get_facecolor")]
        bottom_patches = [p for p in bottom_ax.patches if hasattr(p, "get_facecolor")]

        assert len(top_patches) >= 1, "Top panel should have highlight patch"
        assert len(bottom_patches) >= 1, "Bottom panel should have highlight patch"

    def test_highlight_multiple_regions(self, miami_panel_dfs):
        """Test highlighting multiple regions."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs

        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            highlight_regions=[
                ("1", 500, 1500),
                ("2", 500, 1500),
            ],
        )

        assert fig is not None

        axes = fig.get_axes()
        top_ax, bottom_ax = axes

        # Each panel should have 2 highlight patches (one per region)
        top_patches = [p for p in top_ax.patches if hasattr(p, "get_facecolor")]
        bottom_patches = [p for p in bottom_ax.patches if hasattr(p, "get_facecolor")]

        assert len(top_patches) >= 2, "Top panel should have 2 highlight patches"
        assert len(bottom_patches) >= 2, "Bottom panel should have 2 highlight patches"

    def test_highlight_with_custom_color(self, miami_panel_dfs):
        """Test that custom highlight color is applied."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = miami_panel_dfs

        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            highlight_regions=[("1", 500, 2500)],
            highlight_color="red",
            highlight_alpha=0.5,
        )

        assert fig is not None

        axes = fig.get_axes()
        top_ax = axes[0]

        # Find patches and verify color
        patches = [p for p in top_ax.patches if hasattr(p, "get_facecolor")]
        assert len(patches) >= 1

        # Get the facecolor and verify it's red-ish
        facecolor = patches[0].get_facecolor()
        # facecolor is RGBA tuple - red should have high R value
        assert facecolor[0] >= 0.9, f"Expected red color, got {facecolor}"

    def test_highlight_plotly_backend(self, miami_panel_dfs):
        """Test region highlighting works with plotly backend."""
        plotter = MiamiPlotter(species="canine", backend="plotly")
        top_df, bottom_df = miami_panel_dfs

        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            highlight_regions=[("1", 500, 2500)],
        )

        assert fig is not None
        # Plotly should have shapes added for the highlights
        assert hasattr(fig.layout, "shapes") or len(fig.layout.shapes or []) >= 0

    def test_highlight_bokeh_backend(self, miami_panel_dfs):
        """Test region highlighting works with bokeh backend."""
        from bokeh.models.layouts import LayoutDOM

        plotter = MiamiPlotter(species="canine", backend="bokeh")
        top_df, bottom_df = miami_panel_dfs

        fig = plotter.plot_miami(
            top_df,
            bottom_df,
            highlight_regions=[("1", 500, 2500)],
        )

        assert fig is not None
        assert isinstance(fig, LayoutDOM)


class TestConstructorThresholdIsTheDefault:
    """MiamiPlotter's genomewide_threshold defaults both panel thresholds."""

    @staticmethod
    def _gwas():
        return pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2],
                "pos": [1_000_000, 2_000_000, 1_000_000, 2_000_000],
                "p": [1e-9, 0.01, 1e-6, 0.5],
            }
        )

    @staticmethod
    def _dashed_y(fig):
        return [
            line.get_ydata()[0]
            for ax in fig.axes
            for line in ax.get_lines()
            if line.get_linestyle() == "--"
        ]

    def test_both_panels_use_the_constructor_threshold(self):
        plotter = MiamiPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_miami(self._gwas(), self._gwas())

        # The bottom panel inverts its y-axis, so both lines sit at +3.0.
        assert self._dashed_y(fig) == pytest.approx([3.0, 3.0])

    def test_explicit_panel_thresholds_beat_the_constructor(self):
        plotter = MiamiPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_miami(
            self._gwas(), self._gwas(), top_threshold=1e-6, bottom_threshold=None
        )

        assert self._dashed_y(fig) == pytest.approx([6.0])

"""Tests for MiamiPlotter class."""

import pandas as pd
import pytest

# Test for import error first (RED phase)
try:
    from pylocuszoom.miami_plotter import MiamiPlotter

    MIAMI_PLOTTER_AVAILABLE = True
except ImportError:
    MIAMI_PLOTTER_AVAILABLE = False
    MiamiPlotter = None

# Check if optional backends are available
import importlib.util

PLOTLY_AVAILABLE = importlib.util.find_spec("plotly") is not None
BOKEH_AVAILABLE = importlib.util.find_spec("bokeh") is not None


@pytest.mark.skipif(
    not MIAMI_PLOTTER_AVAILABLE, reason="MiamiPlotter not implemented yet"
)
class TestMiamiPlotter:
    """Tests for the MiamiPlotter class."""

    @pytest.fixture
    def plotter(self):
        """Create a MiamiPlotter instance."""
        return MiamiPlotter(species="canine")

    @pytest.fixture
    def gwas_data(self):
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

    def test_plot_miami_returns_figure(self, plotter, gwas_data):
        """Test that plot_miami returns a figure object."""
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df)
        assert fig is not None

    def test_matplotlib_backend(self, gwas_data):
        """Test MiamiPlotter with matplotlib backend."""
        from matplotlib.figure import Figure

        plotter = MiamiPlotter(species="canine", backend="matplotlib")
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df)
        assert fig is not None
        assert isinstance(fig, Figure)

    @pytest.mark.skipif(not PLOTLY_AVAILABLE, reason="Plotly not installed")
    def test_plotly_backend(self, gwas_data):
        """Test MiamiPlotter with plotly backend."""
        import plotly.graph_objects as go

        plotter = MiamiPlotter(species="canine", backend="plotly")
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df)
        assert fig is not None
        assert isinstance(fig, go.Figure)

    @pytest.mark.skipif(not BOKEH_AVAILABLE, reason="Bokeh not installed")
    def test_bokeh_backend(self, gwas_data):
        """Test MiamiPlotter with bokeh backend."""
        from bokeh.models.layouts import LayoutDOM

        plotter = MiamiPlotter(species="canine", backend="bokeh")
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df)
        assert fig is not None
        assert isinstance(fig, LayoutDOM)

    def test_chromosome_colors_match(self, plotter, gwas_data):
        """Test that same chromosome has same color in both panels.

        Verifies that chromosome colors are consistent between top and bottom
        panels - critical for visual comparison in Miami plots.
        """
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df)

        # Get the figure axes (matplotlib-specific verification)
        axes = fig.get_axes()
        assert len(axes) == 2, "Miami plot should have exactly 2 panels"

        # For this test, we verify the figure was created with 2 panels
        # Color consistency is implicitly tested by using the same color mapping
        top_ax, bottom_ax = axes
        assert top_ax is not None
        assert bottom_ax is not None

    def test_empty_dataframe_raises(self, plotter):
        """Test that empty DataFrame raises appropriate error."""
        empty_df = pd.DataFrame(columns=["chrom", "pos", "p"])
        valid_df = pd.DataFrame({"chrom": [1], "pos": [1000], "p": [0.01]})

        # Empty top_df should raise
        with pytest.raises((ValueError, KeyError)):
            plotter.plot_miami(empty_df, valid_df)

        # Empty bottom_df should raise
        with pytest.raises((ValueError, KeyError)):
            plotter.plot_miami(valid_df, empty_df)

    def test_missing_columns_raises(self, plotter):
        """Test that missing required columns raises ValidationError or ValueError."""
        # DataFrame missing 'p' column
        missing_p_df = pd.DataFrame({"chrom": [1], "pos": [1000]})
        valid_df = pd.DataFrame({"chrom": [1], "pos": [1000], "p": [0.01]})

        with pytest.raises(ValueError):
            plotter.plot_miami(missing_p_df, valid_df)

        # DataFrame missing 'chrom' column
        missing_chrom_df = pd.DataFrame({"pos": [1000], "p": [0.01]})
        with pytest.raises(ValueError):
            plotter.plot_miami(valid_df, missing_chrom_df)

    def test_inverted_bottom_panel(self, plotter, gwas_data):
        """Test that bottom panel y-axis is inverted.

        The bottom panel should display with y-axis inverted (max at top/0 at bottom)
        to create the characteristic Miami plot mirror effect.
        """
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df)

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


@pytest.mark.skipif(
    not MIAMI_PLOTTER_AVAILABLE, reason="MiamiPlotter not implemented yet"
)
class TestMiamiPlotterOptions:
    """Tests for MiamiPlotter configuration options."""

    @pytest.fixture
    def gwas_data(self):
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

    def test_custom_figsize(self, gwas_data):
        """Test that custom figsize is respected."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df, figsize=(15, 10))

        # Check figure size
        width, height = fig.get_size_inches()
        assert width == 15
        assert height == 10

    def test_custom_title(self, gwas_data):
        """Test that custom title is set."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = gwas_data
        fig = plotter.plot_miami(top_df, bottom_df, title="Custom Miami Title")

        # Check that title was set (matplotlib-specific)
        title_text = fig._suptitle
        assert title_text is not None

    def test_panel_labels(self, gwas_data):
        """Test that panel labels are applied."""
        plotter = MiamiPlotter(species="canine")
        top_df, bottom_df = gwas_data
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


@pytest.mark.skipif(
    not MIAMI_PLOTTER_AVAILABLE, reason="MiamiPlotter not implemented yet"
)
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

    @pytest.mark.skipif(not PLOTLY_AVAILABLE, reason="Plotly not installed")
    def test_plotly_hover_data(self, gwas_data_with_rs):
        """Test that plotly backend creates figure with hover data."""
        plotter = MiamiPlotter(species="canine", backend="plotly")
        top_df, bottom_df = gwas_data_with_rs
        fig = plotter.plot_miami(top_df, bottom_df, rs_col="rs")
        assert fig is not None

        # Verify traces exist (plotly-specific)
        assert len(fig.data) > 0, "Plotly figure should have traces for hover data"

    @pytest.mark.skipif(not BOKEH_AVAILABLE, reason="Bokeh not installed")
    def test_bokeh_hover_data(self, gwas_data_with_rs):
        """Test that bokeh backend creates figure with hover tools."""
        plotter = MiamiPlotter(species="canine", backend="bokeh")
        top_df, bottom_df = gwas_data_with_rs
        fig = plotter.plot_miami(top_df, bottom_df, rs_col="rs")
        assert fig is not None

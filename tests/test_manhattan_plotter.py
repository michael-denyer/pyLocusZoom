"""Tests for ManhattanPlotter class."""

import pandas as pd
import pytest

from pylocuszoom.manhattan_plotter import ManhattanPlotter


class TestManhattanPlotter:
    """Tests for the ManhattanPlotter class."""

    @pytest.fixture
    def plotter(self):
        """Create a ManhattanPlotter instance."""
        return ManhattanPlotter(species="canine")

    @pytest.fixture
    def gwas_data(self):
        """Create sample GWAS data."""
        return pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3],
                "pos": [1000, 2000, 1000, 2000, 1000],
                "p": [0.01, 0.001, 0.0001, 0.05, 1e-8],
            }
        )

    def test_plot_manhattan_returns_figure(self, plotter, gwas_data):
        """Test that plot_manhattan returns a figure object."""
        fig = plotter.plot_manhattan(gwas_data)
        assert fig is not None

    def test_plot_qq_returns_figure(self, plotter, gwas_data):
        """Test that plot_qq returns a figure object."""
        fig = plotter.plot_qq(gwas_data)
        assert fig is not None

    def test_plot_manhattan_qq_returns_figure(self, plotter, gwas_data):
        """Test that plot_manhattan_qq returns a figure object."""
        fig = plotter.plot_manhattan_qq(gwas_data)
        assert fig is not None

    def test_plot_manhattan_stacked_returns_figure(self, plotter, gwas_data):
        """Test that plot_manhattan_stacked returns a figure object."""
        fig = plotter.plot_manhattan_stacked([gwas_data, gwas_data])
        assert fig is not None


class TestManhattanPlotterBackends:
    """Tests for ManhattanPlotter backend support."""

    @pytest.fixture
    def gwas_data(self):
        """Create sample GWAS data."""
        return pd.DataFrame(
            {
                "chrom": [1, 1, 2],
                "pos": [1000, 2000, 1000],
                "p": [0.01, 0.001, 0.0001],
            }
        )

    def test_matplotlib_backend(self, gwas_data):
        """Test ManhattanPlotter with matplotlib backend."""
        plotter = ManhattanPlotter(species="canine", backend="matplotlib")
        fig = plotter.plot_manhattan(gwas_data)
        assert fig is not None

    @pytest.mark.skipif(
        True,
        reason="Plotly backend requires plotly package",
    )
    def test_plotly_backend(self, gwas_data):
        """Test ManhattanPlotter with plotly backend."""
        plotter = ManhattanPlotter(species="canine", backend="plotly")
        fig = plotter.plot_manhattan(gwas_data)
        assert fig is not None


class TestCategoricalManhattanIntegerCategories:
    """Test categorical Manhattan with integer category columns."""

    def test_integer_categories_render_all_points(self):
        """Categorical Manhattan should render all points, not silently drop by type mismatch."""
        import matplotlib.pyplot as plt

        df = pd.DataFrame(
            {
                "cat": [1, 1, 2, 2, 3],
                "p": [0.01, 0.05, 0.1, 0.001, 0.5],
            }
        )
        plotter = ManhattanPlotter(species="human")
        fig = plotter.plot_manhattan(df, category_col="cat", p_col="p")
        assert fig is not None

        ax = fig.axes[0]
        total_points = sum(len(coll.get_offsets()) for coll in ax.collections)
        assert total_points == 5, f"Expected 5 points rendered, got {total_points}"
        plt.close(fig)

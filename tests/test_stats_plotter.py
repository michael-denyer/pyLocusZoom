"""Tests for StatsPlotter class."""

import pandas as pd
import pytest

from pylocuszoom.stats_plotter import StatsPlotter


class TestStatsPlotter:
    """Tests for the StatsPlotter class."""

    @pytest.fixture
    def plotter(self):
        """Create a StatsPlotter instance."""
        return StatsPlotter()

    @pytest.fixture
    def phewas_data(self):
        """Create sample PheWAS data."""
        return pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "BMI"],
                "category": ["Anthropometric", "Anthropometric", "Anthropometric"],
                "p_value": [0.01, 0.001, 1e-8],
            }
        )

    @pytest.fixture
    def forest_data(self):
        """Create sample forest plot data."""
        return pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [0.5, 0.3, 0.4],
                "ci_lower": [0.2, 0.1, 0.2],
                "ci_upper": [0.8, 0.5, 0.6],
            }
        )

    def test_plot_phewas_returns_figure(self, plotter, phewas_data):
        """Test that plot_phewas returns a figure object."""
        fig = plotter.plot_phewas(phewas_data, variant_id="rs12345")
        assert fig is not None

    def test_plot_forest_returns_figure(self, plotter, forest_data):
        """Test that plot_forest returns a figure object."""
        fig = plotter.plot_forest(forest_data, variant_id="rs12345")
        assert fig is not None


class TestStatsPlotterBackends:
    """Tests for StatsPlotter backend support."""

    @pytest.fixture
    def phewas_data(self):
        """Create sample PheWAS data."""
        return pd.DataFrame(
            {
                "phenotype": ["Height", "Weight"],
                "category": ["Anthro", "Anthro"],
                "p_value": [0.01, 0.001],
            }
        )

    def test_matplotlib_backend(self, phewas_data):
        """Test StatsPlotter with matplotlib backend."""
        plotter = StatsPlotter(backend="matplotlib")
        fig = plotter.plot_phewas(phewas_data, variant_id="rs12345")
        assert fig is not None

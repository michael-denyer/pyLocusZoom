"""Tests for StatsPlotter class."""

import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import PathCollection
from matplotlib.colors import to_hex
from matplotlib.markers import MarkerStyle

from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.colors import get_phewas_category_color
from pylocuszoom.stats_plotter import StatsPlotter
from tests.conftest import FIGURE_TYPES


def _scatters(fig):
    """Return the scatter collections drawn on the figure's first axes."""
    return [c for c in fig.get_axes()[0].collections if isinstance(c, PathCollection)]


def _marker_name(collection):
    """Name the built-in marker whose path a scatter collection draws."""
    drawn = collection.get_paths()[0].vertices
    for name in ("o", "^", "v", "s"):
        style = MarkerStyle(name)
        expected = style.get_path().transformed(style.get_transform()).vertices
        if drawn.shape == expected.shape and np.allclose(drawn, expected):
            return name
    raise AssertionError(f"unrecognised marker path: {drawn.tolist()}")


def _markers(fig):
    """Pair each scatter's marker shape with the number of points it draws."""
    return [(_marker_name(c), len(c.get_offsets())) for c in _scatters(fig)]


def _colour_counts(fig):
    """Pair each scatter's face colour with the number of points it draws."""
    return [
        (to_hex(c.get_facecolor()[0]), len(c.get_offsets())) for c in _scatters(fig)
    ]


def _dashed_x(fig):
    """Return the x position of every dashed reference line on the figure."""
    return [
        line.get_xdata()[0]
        for ax in fig.axes
        for line in ax.get_lines()
        if line.get_linestyle() == "--"
    ]


class TestStatsPlotter:
    """Tests for the StatsPlotter class."""

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

    def test_plot_phewas_returns_figure(self, stats_plotter, phewas_data):
        """Return a matplotlib figure carrying one point per phenotype."""
        fig = stats_plotter.plot_phewas(phewas_data, variant_id="rs12345")

        assert isinstance(fig, FIGURE_TYPES["matplotlib"])
        assert _markers(fig) == [("o", len(phewas_data))]

    def test_plot_forest_returns_figure(self, stats_plotter, forest_data):
        """Return a matplotlib figure marking every study at its effect size."""
        fig = stats_plotter.plot_forest(forest_data, variant_id="rs12345")

        assert isinstance(fig, FIGURE_TYPES["matplotlib"])
        (collection,) = _scatters(fig)
        assert collection.get_offsets()[:, 0].tolist() == pytest.approx(
            forest_data["effect"].tolist()
        )


class TestPheWASEdgeCases:
    """Tests for PheWAS plot edge cases."""

    def test_phewas_without_category_column(self, stats_plotter):
        """Draw every phenotype in one ungrouped scatter when no category exists."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "BMI"],
                "p_value": [0.01, 0.001, 1e-8],
            }
        )

        fig = stats_plotter.plot_phewas(
            df,
            variant_id="rs12345",
            phenotype_col="phenotype",
            p_col="p_value",
            category_col="nonexistent",
        )

        assert _markers(fig) == [("o", 3)]

    def test_phewas_with_effect_column_positive(self, stats_plotter):
        """Draw only upward triangles when every effect is positive."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight"],
                "category": ["Anthro", "Anthro"],
                "p_value": [0.01, 0.001],
                "beta": [0.5, 0.3],
            }
        )

        fig = stats_plotter.plot_phewas(df, variant_id="rs12345", effect_col="beta")

        assert _markers(fig) == [("^", 2)]

    def test_phewas_with_effect_column_negative(self, stats_plotter):
        """Draw only downward triangles when every effect is negative."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight"],
                "category": ["Anthro", "Anthro"],
                "p_value": [0.01, 0.001],
                "beta": [-0.5, -0.3],
            }
        )

        fig = stats_plotter.plot_phewas(df, variant_id="rs12345", effect_col="beta")

        assert _markers(fig) == [("v", 2)]

    def test_phewas_with_mixed_effects(self, stats_plotter):
        """Split each category into an upward and a downward triangle scatter."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "BMI", "HDL"],
                "category": ["Anthro", "Anthro", "Anthro", "Lipids"],
                "p_value": [0.01, 0.001, 1e-5, 1e-6],
                "beta": [0.5, -0.3, 0.2, -0.4],
            }
        )

        fig = stats_plotter.plot_phewas(df, variant_id="rs12345", effect_col="beta")

        assert _markers(fig) == [("^", 2), ("v", 1), ("v", 1)]

    def test_phewas_with_multiple_categories(self, stats_plotter):
        """Give each category its own scatter in the next palette colour."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "LDL", "HDL", "Glucose"],
                "category": ["Anthro", "Anthro", "Lipids", "Lipids", "Metabolic"],
                "p_value": [0.01, 0.001, 1e-5, 1e-6, 1e-4],
            }
        )

        fig = stats_plotter.plot_phewas(df, variant_id="rs12345")

        assert _colour_counts(fig) == [
            (get_phewas_category_color(0).lower(), 2),
            (get_phewas_category_color(1).lower(), 2),
            (get_phewas_category_color(2).lower(), 1),
        ]

    def test_phewas_with_nan_category(self, stats_plotter):
        """A null category is its own group and takes the next palette colour."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "Unknown"],
                "category": ["Anthro", "Anthro", None],  # NaN category
                "p_value": [0.01, 0.001, 1e-5],
            }
        )

        fig = stats_plotter.plot_phewas(df, variant_id="rs12345")

        assert _colour_counts(fig) == [
            (get_phewas_category_color(0).lower(), 2),
            (get_phewas_category_color(1).lower(), 1),
        ]


class TestForestPlotEdgeCases:
    """Tests for forest plot edge cases."""

    def test_forest_with_weight_column(self, stats_plotter):
        """Scale each marker between 40 and 200 across the weight range."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [0.5, 0.3, 0.4],
                "ci_lower": [0.2, 0.1, 0.2],
                "ci_upper": [0.8, 0.5, 0.6],
                "weight": [10.0, 30.0, 20.0],
            }
        )

        fig = stats_plotter.plot_forest(df, variant_id="rs12345", weight_col="weight")

        (collection,) = _scatters(fig)
        assert collection.get_sizes().tolist() == pytest.approx([40.0, 200.0, 120.0])

    def test_forest_with_equal_weights(self, stats_plotter):
        """Give every marker one shared size when the weights do not vary."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [0.5, 0.3, 0.4],
                "ci_lower": [0.2, 0.1, 0.2],
                "ci_upper": [0.8, 0.5, 0.6],
                "weight": [20.0, 20.0, 20.0],
            }
        )

        fig = stats_plotter.plot_forest(df, variant_id="rs12345", weight_col="weight")

        (collection,) = _scatters(fig)
        assert collection.get_sizes().tolist() == pytest.approx([120.0])

    def test_forest_with_odds_ratio(self, stats_plotter):
        """Put the null line at 1.0 and label the x axis Odds Ratio."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [1.5, 0.8, 1.2],
                "ci_lower": [1.1, 0.5, 0.9],
                "ci_upper": [2.0, 1.2, 1.6],
            }
        )

        fig = stats_plotter.plot_forest(
            df,
            variant_id="rs12345",
            null_value=1.0,
            effect_label="Odds Ratio",
        )

        assert _dashed_x(fig) == pytest.approx([1.0])
        assert fig.get_axes()[0].get_xlabel() == "Odds Ratio"

    def test_forest_with_custom_figsize(self, stats_plotter):
        """Size the figure to the requested width and height."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B"],
                "effect": [0.5, 0.3],
                "ci_lower": [0.2, 0.1],
                "ci_upper": [0.8, 0.5],
            }
        )

        fig = stats_plotter.plot_forest(df, variant_id="rs12345", figsize=(12, 10))

        assert fig.get_size_inches().tolist() == pytest.approx([12.0, 10.0])


class TestStatsPlotterInit:
    """Tests for StatsPlotter initialization."""

    def test_default_initialization(self):
        """StatsPlotter should initialize with default parameters."""
        plotter = StatsPlotter()
        assert plotter.genomewide_threshold == 5e-8

    def test_custom_genomewide_threshold(self):
        """StatsPlotter should accept custom genomewide threshold."""
        plotter = StatsPlotter(genomewide_threshold=1e-5)
        assert plotter.genomewide_threshold == 1e-5

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_accepts_every_builtin_backend(self, backend):
        """Return that backend's own figure type from plot_phewas."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight"],
                "category": ["Anthro", "Anthro"],
                "p_value": [0.01, 0.001],
            }
        )

        fig = StatsPlotter(backend=backend).plot_phewas(df, variant_id="rs12345")

        assert isinstance(fig, FIGURE_TYPES[backend])


class TestConstructorThresholdIsTheDefault:
    """StatsPlotter's genomewide_threshold reaches plot_phewas."""

    @staticmethod
    def _phewas():
        return pd.DataFrame(
            {
                "phenotype": ["a", "b", "c"],
                "p_value": [1e-9, 1e-4, 0.3],
                "category": ["x", "x", "y"],
            }
        )

    def test_plot_phewas_uses_the_constructor_threshold(self):
        plotter = StatsPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_phewas(self._phewas(), variant_id="rs1")

        assert _dashed_x(fig) == pytest.approx([3.0])

    def test_explicit_argument_beats_the_constructor(self):
        plotter = StatsPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_phewas(
            self._phewas(), variant_id="rs1", significance_threshold=1e-6
        )

        assert _dashed_x(fig) == pytest.approx([6.0])

    def test_explicit_none_still_draws_no_line(self):
        plotter = StatsPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_phewas(
            self._phewas(), variant_id="rs1", significance_threshold=None
        )

        assert _dashed_x(fig) == []


class TestPheWASManyCategories:
    """Test PheWAS plot with many categories cycles colors correctly."""

    def test_phewas_15_categories_succeeds(self, stats_plotter):
        """PheWAS with >12 categories should cycle colors without error.

        The PHEWAS_CATEGORY_COLORS palette has 12 colors, so 15 categories
        requires color cycling. This tests that modulo indexing works.
        """
        # Create 15 unique categories
        categories = [f"Category_{i}" for i in range(15)]
        phenotypes = [f"Phenotype_{i}" for i in range(15)]

        phewas_df = pd.DataFrame(
            {
                "phenotype": phenotypes,
                "p_value": [10 ** (-i - 1) for i in range(15)],  # varying p-values
                "category": categories,
            }
        )

        fig = stats_plotter.plot_phewas(
            phewas_df,
            variant_id="rs12345",
            phenotype_col="phenotype",
            p_col="p_value",
            category_col="category",
        )

        ax = fig.axes[0]
        total_points = sum(
            len(collection.get_offsets()) for collection in ax.collections
        )
        assert total_points == 15

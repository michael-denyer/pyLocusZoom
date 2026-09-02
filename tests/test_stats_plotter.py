"""Tests for StatsPlotter class."""

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.stats_plotter import StatsPlotter


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
        """Test that plot_phewas returns a figure object."""
        fig = stats_plotter.plot_phewas(phewas_data, variant_id="rs12345")
        assert fig is not None
        plt.close(fig)

    def test_plot_forest_returns_figure(self, stats_plotter, forest_data):
        """Test that plot_forest returns a figure object."""
        fig = stats_plotter.plot_forest(forest_data, variant_id="rs12345")
        assert fig is not None
        plt.close(fig)


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
        plt.close(fig)


class TestPheWASEdgeCases:
    """Tests for PheWAS plot edge cases."""

    def test_phewas_without_category_column(self, stats_plotter):
        """PheWAS plot should work without category column."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "BMI"],
                "p_value": [0.01, 0.001, 1e-8],
            }
        )
        # Remove category column entirely - test the else branch at line 92-94
        fig = stats_plotter.plot_phewas(
            df,
            variant_id="rs12345",
            phenotype_col="phenotype",
            p_col="p_value",
            category_col="nonexistent",  # Column doesn't exist
        )
        assert fig is not None
        plt.close(fig)

    def test_phewas_with_effect_column_positive(self, stats_plotter):
        """PheWAS plot with positive effects should show upward triangles."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight"],
                "category": ["Anthro", "Anthro"],
                "p_value": [0.01, 0.001],
                "beta": [0.5, 0.3],  # All positive effects
            }
        )
        fig = stats_plotter.plot_phewas(df, variant_id="rs12345", effect_col="beta")
        assert fig is not None
        plt.close(fig)

    def test_phewas_with_effect_column_negative(self, stats_plotter):
        """PheWAS plot with negative effects should show downward triangles."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight"],
                "category": ["Anthro", "Anthro"],
                "p_value": [0.01, 0.001],
                "beta": [-0.5, -0.3],  # All negative effects
            }
        )
        fig = stats_plotter.plot_phewas(df, variant_id="rs12345", effect_col="beta")
        assert fig is not None
        plt.close(fig)

    def test_phewas_with_mixed_effects(self, stats_plotter):
        """PheWAS plot with mixed effects should show both markers."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "BMI", "HDL"],
                "category": ["Anthro", "Anthro", "Anthro", "Lipids"],
                "p_value": [0.01, 0.001, 1e-5, 1e-6],
                "beta": [0.5, -0.3, 0.2, -0.4],  # Mixed effects
            }
        )
        fig = stats_plotter.plot_phewas(df, variant_id="rs12345", effect_col="beta")
        assert fig is not None
        plt.close(fig)

    def test_phewas_with_multiple_categories(self, stats_plotter):
        """PheWAS plot with multiple categories should color correctly."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "LDL", "HDL", "Glucose"],
                "category": ["Anthro", "Anthro", "Lipids", "Lipids", "Metabolic"],
                "p_value": [0.01, 0.001, 1e-5, 1e-6, 1e-4],
            }
        )
        fig = stats_plotter.plot_phewas(df, variant_id="rs12345")
        assert fig is not None
        plt.close(fig)

    def test_phewas_with_nan_category(self, stats_plotter):
        """PheWAS plot should handle NaN category values."""
        df = pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "Unknown"],
                "category": ["Anthro", "Anthro", None],  # NaN category
                "p_value": [0.01, 0.001, 1e-5],
            }
        )
        fig = stats_plotter.plot_phewas(df, variant_id="rs12345")
        assert fig is not None
        plt.close(fig)


class TestForestPlotEdgeCases:
    """Tests for forest plot edge cases."""

    def test_forest_with_weight_column(self, stats_plotter):
        """Forest plot with weight column should size markers."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [0.5, 0.3, 0.4],
                "ci_lower": [0.2, 0.1, 0.2],
                "ci_upper": [0.8, 0.5, 0.6],
                "weight": [10.0, 30.0, 20.0],  # Different weights
            }
        )
        fig = stats_plotter.plot_forest(df, variant_id="rs12345", weight_col="weight")
        assert fig is not None
        plt.close(fig)

    def test_forest_with_equal_weights(self, stats_plotter):
        """Forest plot with equal weights should have uniform marker sizes."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [0.5, 0.3, 0.4],
                "ci_lower": [0.2, 0.1, 0.2],
                "ci_upper": [0.8, 0.5, 0.6],
                "weight": [20.0, 20.0, 20.0],  # Equal weights - triggers line 257-258
            }
        )
        fig = stats_plotter.plot_forest(df, variant_id="rs12345", weight_col="weight")
        assert fig is not None
        plt.close(fig)

    def test_forest_with_odds_ratio(self, stats_plotter):
        """Forest plot with odds ratio should use null_value=1.0."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C"],
                "effect": [1.5, 0.8, 1.2],  # Odds ratios
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
        assert fig is not None
        plt.close(fig)

    def test_forest_with_custom_figsize(self, stats_plotter):
        """Forest plot should accept custom figure size."""
        df = pd.DataFrame(
            {
                "study": ["Study A", "Study B"],
                "effect": [0.5, 0.3],
                "ci_lower": [0.2, 0.1],
                "ci_upper": [0.8, 0.5],
            }
        )
        fig = stats_plotter.plot_forest(df, variant_id="rs12345", figsize=(12, 10))
        assert fig is not None
        plt.close(fig)


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
        """StatsPlotter constructs against any built-in backend."""
        plotter = StatsPlotter(backend=backend)
        assert plotter._backend is not None


class TestBarChartCapabilityGate:
    """Forest plots need error bars; PheWAS plots do not."""

    @staticmethod
    def _forest_df():
        return pd.DataFrame(
            {
                "study": ["A", "B"],
                "beta": [0.1, 0.2],
                "ci_lower": [0.0, 0.1],
                "ci_upper": [0.2, 0.3],
            }
        )

    def test_forest_rejects_backend_without_error_bars(self):
        from pylocuszoom._stats_renderer import StatsRenderer

        class BarelyABackend:
            pass

        renderer = StatsRenderer(BarelyABackend())

        with pytest.raises(TypeError, match="does not support bar charts"):
            renderer.render_forest(
                self._forest_df(),
                variant_id="rs1",
                study_col="study",
                effect_col="beta",
                ci_lower_col="ci_lower",
                ci_upper_col="ci_upper",
                weight_col=None,
                null_value=0.0,
                effect_label="Beta",
                figsize=(6.0, 4.0),
            )

    def _render(self, backend):
        from pylocuszoom._stats_renderer import StatsRenderer

        StatsRenderer(backend).render_forest(
            self._forest_df(),
            variant_id="rs1",
            study_col="study",
            effect_col="beta",
            ci_lower_col="ci_lower",
            ci_upper_col="ci_upper",
            weight_col=None,
            null_value=0.0,
            effect_label="Beta",
            figsize=(6.0, 4.0),
        )

    def test_error_names_the_backend_and_the_protocol(self):
        class BarelyABackend:
            pass

        with pytest.raises(TypeError) as excinfo:
            self._render(BarelyABackend())

        assert "BarelyABackend" in str(excinfo.value)
        assert "SupportsBarCharts" in str(excinfo.value)

    def test_error_names_both_required_methods(self):
        """errorbar_h alone does not satisfy the protocol, so say what is missing."""

        class ErrorBarOnlyBackend:
            def errorbar_h(self, *args, **kwargs): ...

        with pytest.raises(TypeError) as excinfo:
            self._render(ErrorBarOnlyBackend())

        message = str(excinfo.value)
        assert "hbar" in message
        assert "errorbar_h" in message


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

    @staticmethod
    def _dashed_x(fig):
        return [
            line.get_xdata()[0]
            for ax in fig.axes
            for line in ax.get_lines()
            if line.get_linestyle() == "--"
        ]

    def test_plot_phewas_uses_the_constructor_threshold(self):
        plotter = StatsPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_phewas(self._phewas(), variant_id="rs1")

        assert self._dashed_x(fig) == pytest.approx([3.0])

    def test_explicit_argument_beats_the_constructor(self):
        plotter = StatsPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_phewas(
            self._phewas(), variant_id="rs1", significance_threshold=1e-6
        )

        assert self._dashed_x(fig) == pytest.approx([6.0])

    def test_explicit_none_still_draws_no_line(self):
        plotter = StatsPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_phewas(
            self._phewas(), variant_id="rs1", significance_threshold=None
        )

        assert self._dashed_x(fig) == []

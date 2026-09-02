"""Tests for ManhattanPlotter class."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.manhattan_plotter import ManhattanPlotter


class TestManhattanPlotter:
    """Tests for the ManhattanPlotter class."""

    @pytest.fixture
    def canine_manhattan_plotter(self):
        """Create a ManhattanPlotter instance."""
        return ManhattanPlotter(species="canine")

    @pytest.fixture
    def manhattan_chrom_df(self):
        """Create sample GWAS data."""
        return pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2, 3],
                "pos": [1000, 2000, 1000, 2000, 1000],
                "p": [0.01, 0.001, 0.0001, 0.05, 1e-8],
            }
        )

    def test_plot_manhattan_returns_figure(
        self, canine_manhattan_plotter, manhattan_chrom_df
    ):
        """Test that plot_manhattan returns a figure object."""
        fig = canine_manhattan_plotter.plot_manhattan(manhattan_chrom_df)
        assert fig is not None

    def test_plot_qq_returns_figure(self, canine_manhattan_plotter, manhattan_chrom_df):
        """Test that plot_qq returns a figure object."""
        fig = canine_manhattan_plotter.plot_qq(manhattan_chrom_df)
        assert fig is not None

    def test_plot_manhattan_qq_returns_figure(
        self, canine_manhattan_plotter, manhattan_chrom_df
    ):
        """Test that plot_manhattan_qq returns a figure object."""
        fig = canine_manhattan_plotter.plot_manhattan_qq(manhattan_chrom_df)
        assert fig is not None

    def test_plot_manhattan_stacked_returns_figure(
        self, canine_manhattan_plotter, manhattan_chrom_df
    ):
        """Test that plot_manhattan_stacked returns a figure object."""
        fig = canine_manhattan_plotter.plot_manhattan_stacked(
            [manhattan_chrom_df, manhattan_chrom_df]
        )
        assert fig is not None


class TestManhattanPlotterBackends:
    """Tests for ManhattanPlotter backend support."""

    @pytest.fixture
    def manhattan_chrom_df(self):
        """Create sample GWAS data."""
        return pd.DataFrame(
            {
                "chrom": [1, 1, 2],
                "pos": [1000, 2000, 1000],
                "p": [0.01, 0.001, 0.0001],
            }
        )

    def test_matplotlib_backend(self, manhattan_chrom_df):
        """Test ManhattanPlotter with matplotlib backend."""
        plotter = ManhattanPlotter(species="canine", backend="matplotlib")
        fig = plotter.plot_manhattan(manhattan_chrom_df)
        assert fig is not None

    def test_plotly_backend(self, manhattan_chrom_df):
        """Test ManhattanPlotter with plotly backend."""
        plotter = ManhattanPlotter(species="canine", backend="plotly")
        fig = plotter.plot_manhattan(manhattan_chrom_df)
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


class TestConstructorThresholdIsTheDefault:
    """The constructor's genomewide_threshold reaches every plot method."""

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

    def test_plot_manhattan_uses_the_constructor_threshold(self):
        fig = ManhattanPlotter(genomewide_threshold=1e-3).plot_manhattan(self._gwas())

        assert self._dashed_y(fig) == pytest.approx([3.0])

    def test_plot_manhattan_stacked_uses_the_constructor_threshold(self):
        plotter = ManhattanPlotter(genomewide_threshold=1e-3)
        fig = plotter.plot_manhattan_stacked([self._gwas(), self._gwas()])

        assert self._dashed_y(fig) == pytest.approx([3.0, 3.0])

    def test_plot_manhattan_qq_uses_the_constructor_threshold(self):
        fig = ManhattanPlotter(genomewide_threshold=1e-3).plot_manhattan_qq(
            self._gwas()
        )

        assert any(y == pytest.approx(3.0) for y in self._dashed_y(fig))

    def test_categorical_manhattan_uses_the_constructor_threshold(self):
        df = self._gwas().assign(category=["a", "a", "b", "b"])
        fig = ManhattanPlotter(genomewide_threshold=1e-3).plot_manhattan(
            df, category_col="category"
        )

        assert self._dashed_y(fig) == pytest.approx([3.0])

    def test_explicit_argument_beats_the_constructor(self):
        fig = ManhattanPlotter(genomewide_threshold=1e-3).plot_manhattan(
            self._gwas(), significance_threshold=1e-6
        )

        assert self._dashed_y(fig) == pytest.approx([6.0])

    def test_explicit_none_still_draws_no_line(self):
        fig = ManhattanPlotter(genomewide_threshold=1e-3).plot_manhattan(
            self._gwas(), significance_threshold=None
        )

        assert self._dashed_y(fig) == []

    def test_default_construction_keeps_the_5e_8_line(self):
        fig = ManhattanPlotter().plot_manhattan(self._gwas())

        assert self._dashed_y(fig) == pytest.approx([-np.log10(5e-8)])


class TestManhattanSingleChromosome:
    """Test Manhattan plot with single chromosome data."""

    def test_manhattan_single_chromosome_succeeds(self, manhattan_plotter):
        """Manhattan plot with only one chromosome should work."""
        df = pd.DataFrame(
            {
                "chrom": ["1"] * 10,
                "pos": list(range(1000000, 11000000, 1000000)),
                "p": [1e-8, 0.05, 0.01, 1e-6, 0.1, 0.001, 1e-10, 0.5, 0.005, 1e-3],
            }
        )

        fig = manhattan_plotter.plot_manhattan(
            df,
            chrom_col="chrom",
            pos_col="pos",
            p_col="p",
        )

        assert fig is not None
        # Verify x-axis has chromosome label
        ax = fig.axes[0]
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        assert "1" in tick_labels
        plt.close(fig)


class TestManhattanStackedValidation:
    """Test validation in stacked Manhattan plots."""

    @pytest.fixture
    def sample_df(self):
        """Sample GWAS DataFrame."""
        return pd.DataFrame(
            {
                "chrom": ["1", "1", "2", "2"],
                "pos": [1000000, 2000000, 1500000, 3000000],
                "p": [1e-8, 0.05, 0.01, 1e-6],
            }
        )

    def test_manhattan_stacked_empty_list_raises(self, manhattan_plotter):
        """Empty list of GWAS DataFrames should raise ValueError."""
        with pytest.raises(ValueError, match="At least one GWAS DataFrame"):
            manhattan_plotter.plot_manhattan_stacked([])

    def test_manhattan_stacked_mismatched_panel_labels_raises(
        self, manhattan_plotter, sample_df
    ):
        """Mismatched panel_labels length should raise ValueError."""
        gwas_dfs = [sample_df, sample_df.copy()]
        panel_labels = ["Only One"]

        with pytest.raises(ValueError, match="panel_labels"):
            manhattan_plotter.plot_manhattan_stacked(
                gwas_dfs,
                panel_labels=panel_labels,
            )


class TestManhattanQQStackedValidation:
    """Test validation in stacked Manhattan+QQ plots."""

    def test_manhattan_qq_stacked_empty_list_raises(self, manhattan_plotter):
        """Empty list of GWAS DataFrames should raise ValueError."""
        with pytest.raises(ValueError, match="At least one GWAS DataFrame"):
            manhattan_plotter.plot_manhattan_qq_stacked([])


class TestEmptyManhattanInput:
    """Empty input has no axis limits to compute."""

    def test_manhattan_with_empty_df_raises(self):
        """Manhattan plot with empty DataFrame raises ValueError.

        Empty DataFrames have no data to compute axis limits,
        resulting in NaN limits which matplotlib cannot handle.
        """
        plotter = ManhattanPlotter(species="human")
        empty_df = pd.DataFrame(columns=["chrom", "pos", "p"])

        with pytest.raises(ValueError, match="(NaN|Inf|cannot)"):
            plotter.plot_manhattan(
                empty_df,
                chrom_col="chrom",
                pos_col="pos",
                p_col="p",
            )

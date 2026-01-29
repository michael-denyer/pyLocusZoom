"""Tests for Manhattan and QQ plot methods in LocusZoomPlotter."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.plotter import LocusZoomPlotter


class TestPlotManhattan:
    """Tests for plot_manhattan method."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS DataFrame for testing."""
        np.random.seed(42)
        n_variants = 100
        return pd.DataFrame(
            {
                "chrom": np.repeat([1, 2, 3], [40, 30, 30]),
                "pos": np.concatenate(
                    [
                        np.sort(np.random.randint(1e6, 1e8, 40)),
                        np.sort(np.random.randint(1e6, 1e8, 30)),
                        np.sort(np.random.randint(1e6, 1e8, 30)),
                    ]
                ),
                "p": np.random.uniform(1e-10, 1, n_variants),
            }
        )

    @pytest.fixture
    def plotter(self):
        """Create a plotter instance."""
        return LocusZoomPlotter(species="human")

    def test_plot_manhattan_returns_figure(self, plotter, sample_gwas_df):
        """plot_manhattan should return a matplotlib figure."""
        fig = plotter.plot_manhattan(sample_gwas_df)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_with_custom_columns(self, plotter):
        """plot_manhattan should work with custom column names."""
        df = pd.DataFrame(
            {
                "chromosome": [1, 1, 2],
                "position": [1e6, 2e6, 1e6],
                "pvalue": [1e-8, 0.01, 0.5],
            }
        )
        fig = plotter.plot_manhattan(
            df, chrom_col="chromosome", pos_col="position", p_col="pvalue"
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_with_species_order(self, sample_gwas_df):
        """plot_manhattan should use species-specific chromosome order."""
        plotter = LocusZoomPlotter(species="canine")
        fig = plotter.plot_manhattan(sample_gwas_df)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_with_custom_order(self, plotter, sample_gwas_df):
        """plot_manhattan should accept custom chromosome order."""
        fig = plotter.plot_manhattan(sample_gwas_df, custom_chrom_order=["3", "2", "1"])
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_shows_significance_line(self, plotter, sample_gwas_df):
        """plot_manhattan should show genome-wide significance line by default."""
        fig = plotter.plot_manhattan(sample_gwas_df)
        # Check that a horizontal line exists (either hline or annotation)
        ax = fig.get_axes()[0]
        lines = ax.get_lines()
        # At least one line should be the significance threshold
        assert len(lines) >= 1
        plt.close(fig)

    def test_plot_manhattan_custom_threshold(self, plotter, sample_gwas_df):
        """plot_manhattan should accept custom significance threshold."""
        fig = plotter.plot_manhattan(sample_gwas_df, significance_threshold=1e-5)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_no_threshold(self, plotter, sample_gwas_df):
        """plot_manhattan should allow disabling significance line."""
        fig = plotter.plot_manhattan(sample_gwas_df, significance_threshold=None)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_with_figsize(self, plotter, sample_gwas_df):
        """plot_manhattan should accept figsize parameter."""
        fig = plotter.plot_manhattan(sample_gwas_df, figsize=(12, 4))
        assert fig.get_size_inches()[0] == pytest.approx(12, rel=0.1)
        assert fig.get_size_inches()[1] == pytest.approx(4, rel=0.1)
        plt.close(fig)

    def test_plot_manhattan_with_title(self, plotter, sample_gwas_df):
        """plot_manhattan should accept title parameter."""
        fig = plotter.plot_manhattan(sample_gwas_df, title="Test Manhattan")
        ax = fig.get_axes()[0]
        assert "Test Manhattan" in ax.get_title()
        plt.close(fig)

    def test_plot_manhattan_validates_columns(self, plotter):
        """plot_manhattan should raise on missing columns."""
        df = pd.DataFrame({"wrong": [1], "columns": [2]})
        with pytest.raises(ValueError, match="not found"):
            plotter.plot_manhattan(df)

    def test_plot_manhattan_handles_empty_df(self, plotter):
        """plot_manhattan should raise on empty DataFrame."""
        df = pd.DataFrame({"chrom": [], "pos": [], "p": []})
        # Empty DF causes axis limits to be NaN/Inf
        with pytest.raises((ValueError, Exception)):
            plotter.plot_manhattan(df)

    def test_plot_manhattan_plotly_backend(self, sample_gwas_df):
        """plot_manhattan should work with plotly backend."""
        pytest.importorskip("plotly")
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(species="human", backend="plotly")
        fig = plotter.plot_manhattan(sample_gwas_df)
        assert isinstance(fig, go.Figure)

    def test_plot_manhattan_bokeh_backend(self, sample_gwas_df):
        """plot_manhattan should work with bokeh backend."""
        pytest.importorskip("bokeh")

        plotter = LocusZoomPlotter(species="human", backend="bokeh")
        fig = plotter.plot_manhattan(sample_gwas_df)
        # Bokeh returns a column layout or figure
        assert fig is not None


class TestPlotQQ:
    """Tests for plot_qq method."""

    @pytest.fixture
    def sample_pvalues_df(self):
        """Sample DataFrame with p-values for QQ plot."""
        np.random.seed(42)
        return pd.DataFrame({"p": np.random.uniform(0, 1, 1000)})

    @pytest.fixture
    def plotter(self):
        """Create a plotter instance."""
        return LocusZoomPlotter()

    def test_plot_qq_returns_figure(self, plotter, sample_pvalues_df):
        """plot_qq should return a matplotlib figure."""
        fig = plotter.plot_qq(sample_pvalues_df)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_qq_with_custom_column(self, plotter):
        """plot_qq should work with custom p-value column name."""
        df = pd.DataFrame({"pvalue": np.random.uniform(0, 1, 100)})
        fig = plotter.plot_qq(df, p_col="pvalue")
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_qq_shows_confidence_band(self, plotter, sample_pvalues_df):
        """plot_qq should show confidence band by default."""
        fig = plotter.plot_qq(sample_pvalues_df)
        ax = fig.get_axes()[0]
        # Should have at least 2 artists (points + confidence band fill)
        assert len(ax.collections) >= 1 or len(ax.patches) >= 1
        plt.close(fig)

    def test_plot_qq_no_confidence_band(self, plotter, sample_pvalues_df):
        """plot_qq should allow disabling confidence band."""
        fig = plotter.plot_qq(sample_pvalues_df, show_confidence_band=False)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_qq_shows_diagonal(self, plotter, sample_pvalues_df):
        """plot_qq should show y=x diagonal line."""
        fig = plotter.plot_qq(sample_pvalues_df)
        ax = fig.get_axes()[0]
        lines = ax.get_lines()
        # Should have diagonal line
        assert len(lines) >= 1
        plt.close(fig)

    def test_plot_qq_shows_lambda(self, plotter, sample_pvalues_df):
        """plot_qq should show lambda in title or annotation."""
        fig = plotter.plot_qq(sample_pvalues_df, show_lambda=True)
        ax = fig.get_axes()[0]
        title = ax.get_title()
        # Lambda should be in title or annotation
        assert "λ" in title or "lambda" in title.lower() or len(ax.texts) > 0
        plt.close(fig)

    def test_plot_qq_with_figsize(self, plotter, sample_pvalues_df):
        """plot_qq should accept figsize parameter."""
        fig = plotter.plot_qq(sample_pvalues_df, figsize=(6, 6))
        assert fig.get_size_inches()[0] == pytest.approx(6, rel=0.1)
        assert fig.get_size_inches()[1] == pytest.approx(6, rel=0.1)
        plt.close(fig)

    def test_plot_qq_with_title(self, plotter, sample_pvalues_df):
        """plot_qq should accept title parameter."""
        fig = plotter.plot_qq(sample_pvalues_df, title="Test QQ Plot")
        ax = fig.get_axes()[0]
        assert "Test QQ" in ax.get_title()
        plt.close(fig)

    def test_plot_qq_validates_columns(self, plotter):
        """plot_qq should raise on missing p-value column."""
        df = pd.DataFrame({"wrong": [1, 2, 3]})
        with pytest.raises(ValueError, match="not found"):
            plotter.plot_qq(df)

    def test_plot_qq_handles_all_nan(self, plotter):
        """plot_qq should raise on all NaN p-values."""
        df = pd.DataFrame({"p": [np.nan, np.nan, np.nan]})
        with pytest.raises(ValueError, match="No valid"):
            plotter.plot_qq(df)

    def test_plot_qq_plotly_backend(self, sample_pvalues_df):
        """plot_qq should work with plotly backend."""
        pytest.importorskip("plotly")
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(backend="plotly")
        fig = plotter.plot_qq(sample_pvalues_df)
        assert isinstance(fig, go.Figure)

    def test_plot_qq_bokeh_backend(self, sample_pvalues_df):
        """plot_qq should work with bokeh backend."""
        pytest.importorskip("bokeh")

        plotter = LocusZoomPlotter(backend="bokeh")
        fig = plotter.plot_qq(sample_pvalues_df)
        assert fig is not None


class TestPlotManhattanCategorical:
    """Tests for categorical Manhattan plots (PheWAS-style)."""

    @pytest.fixture
    def sample_phewas_df(self):
        """Sample PheWAS-style DataFrame."""
        return pd.DataFrame(
            {
                "category": ["cardio", "cardio", "neuro", "neuro", "immuno"],
                "phenotype": ["BP", "HR", "AD", "PD", "RA"],
                "p": [1e-10, 0.01, 1e-6, 0.5, 1e-4],
            }
        )

    @pytest.fixture
    def plotter(self):
        """Create a plotter instance."""
        return LocusZoomPlotter()

    def test_plot_manhattan_categorical(self, plotter, sample_phewas_df):
        """plot_manhattan should support categorical x-axis."""
        fig = plotter.plot_manhattan(
            sample_phewas_df,
            category_col="category",
            p_col="p",
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_categorical_custom_order(self, plotter, sample_phewas_df):
        """plot_manhattan should accept custom category order."""
        fig = plotter.plot_manhattan(
            sample_phewas_df,
            category_col="category",
            p_col="p",
            category_order=["neuro", "cardio", "immuno"],
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestPlotManhattanStacked:
    """Tests for stacked Manhattan plots."""

    @pytest.fixture
    def sample_gwas_dfs(self):
        """Multiple sample GWAS DataFrames for stacked testing."""
        np.random.seed(42)
        dfs = []
        for i in range(3):
            n_variants = 50
            dfs.append(
                pd.DataFrame(
                    {
                        "chrom": np.repeat([1, 2], [25, 25]),
                        "pos": np.concatenate(
                            [
                                np.sort(np.random.randint(1e6, 1e8, 25)),
                                np.sort(np.random.randint(1e6, 1e8, 25)),
                            ]
                        ),
                        "p": np.random.uniform(1e-10, 1, n_variants),
                    }
                )
            )
        return dfs

    @pytest.fixture
    def plotter(self):
        """Create a plotter instance."""
        return LocusZoomPlotter(species="human")

    def test_plot_manhattan_stacked_returns_figure(self, plotter, sample_gwas_dfs):
        """plot_manhattan_stacked should return a matplotlib figure."""
        fig = plotter.plot_manhattan_stacked(sample_gwas_dfs)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_stacked_creates_multiple_panels(
        self, plotter, sample_gwas_dfs
    ):
        """plot_manhattan_stacked should create one panel per DataFrame."""
        fig = plotter.plot_manhattan_stacked(sample_gwas_dfs)
        axes = fig.get_axes()
        # Should have 3 panels
        assert len(axes) == 3
        plt.close(fig)

    def test_plot_manhattan_stacked_with_panel_labels(self, plotter, sample_gwas_dfs):
        """plot_manhattan_stacked should show panel labels when provided."""
        labels = ["Study A", "Study B", "Study C"]
        fig = plotter.plot_manhattan_stacked(sample_gwas_dfs, panel_labels=labels)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_stacked_validates_label_count(
        self, plotter, sample_gwas_dfs
    ):
        """plot_manhattan_stacked should raise if panel_labels length mismatch."""
        with pytest.raises(ValueError, match="length"):
            plotter.plot_manhattan_stacked(sample_gwas_dfs, panel_labels=["A", "B"])

    def test_plot_manhattan_stacked_with_figsize(self, plotter, sample_gwas_dfs):
        """plot_manhattan_stacked should accept figsize parameter."""
        fig = plotter.plot_manhattan_stacked(sample_gwas_dfs, figsize=(14, 10))
        assert fig.get_size_inches()[0] == pytest.approx(14, rel=0.1)
        plt.close(fig)

    def test_plot_manhattan_stacked_single_df(self, plotter, sample_gwas_dfs):
        """plot_manhattan_stacked should work with single DataFrame."""
        fig = plotter.plot_manhattan_stacked([sample_gwas_dfs[0]])
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_stacked_empty_list_raises(self, plotter):
        """plot_manhattan_stacked should raise on empty list."""
        with pytest.raises(ValueError, match="At least one"):
            plotter.plot_manhattan_stacked([])

    def test_plot_manhattan_stacked_plotly_backend(self, sample_gwas_dfs):
        """plot_manhattan_stacked should work with plotly backend."""
        pytest.importorskip("plotly")
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(species="human", backend="plotly")
        fig = plotter.plot_manhattan_stacked(sample_gwas_dfs)
        assert isinstance(fig, go.Figure)

    def test_plot_manhattan_stacked_bokeh_backend(self, sample_gwas_dfs):
        """plot_manhattan_stacked should work with bokeh backend."""
        pytest.importorskip("bokeh")

        plotter = LocusZoomPlotter(species="human", backend="bokeh")
        fig = plotter.plot_manhattan_stacked(sample_gwas_dfs)
        assert fig is not None


class TestPlotManhattanQQSideBySide:
    """Tests for side-by-side Manhattan and QQ plots."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS DataFrame for testing."""
        np.random.seed(42)
        n_variants = 100
        return pd.DataFrame(
            {
                "chrom": np.repeat([1, 2, 3], [40, 30, 30]),
                "pos": np.concatenate(
                    [
                        np.sort(np.random.randint(1e6, 1e8, 40)),
                        np.sort(np.random.randint(1e6, 1e8, 30)),
                        np.sort(np.random.randint(1e6, 1e8, 30)),
                    ]
                ),
                "p": np.random.uniform(1e-10, 1, n_variants),
            }
        )

    @pytest.fixture
    def plotter(self):
        """Create a plotter instance."""
        return LocusZoomPlotter(species="human")

    def test_plot_manhattan_qq_returns_figure(self, plotter, sample_gwas_df):
        """plot_manhattan_qq should return a matplotlib figure."""
        fig = plotter.plot_manhattan_qq(sample_gwas_df)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_qq_creates_two_panels(self, plotter, sample_gwas_df):
        """plot_manhattan_qq should create two side-by-side panels."""
        fig = plotter.plot_manhattan_qq(sample_gwas_df)
        axes = fig.get_axes()
        # Should have 2 panels (Manhattan + QQ)
        assert len(axes) == 2
        plt.close(fig)

    def test_plot_manhattan_qq_with_title(self, plotter, sample_gwas_df):
        """plot_manhattan_qq should accept title parameter."""
        fig = plotter.plot_manhattan_qq(sample_gwas_df, title="Combined Plot")
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_qq_with_figsize(self, plotter, sample_gwas_df):
        """plot_manhattan_qq should accept figsize parameter."""
        fig = plotter.plot_manhattan_qq(sample_gwas_df, figsize=(16, 5))
        assert fig.get_size_inches()[0] == pytest.approx(16, rel=0.1)
        plt.close(fig)

    def test_plot_manhattan_qq_custom_columns(self, plotter):
        """plot_manhattan_qq should work with custom column names."""
        df = pd.DataFrame(
            {
                "chromosome": [1, 1, 2],
                "position": [1e6, 2e6, 1e6],
                "pvalue": [1e-8, 0.01, 0.5],
            }
        )
        fig = plotter.plot_manhattan_qq(
            df, chrom_col="chromosome", pos_col="position", p_col="pvalue"
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_manhattan_qq_plotly_backend(self, sample_gwas_df):
        """plot_manhattan_qq should work with plotly backend."""
        pytest.importorskip("plotly")
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(species="human", backend="plotly")
        fig = plotter.plot_manhattan_qq(sample_gwas_df)
        assert isinstance(fig, go.Figure)

    def test_plot_manhattan_qq_bokeh_backend(self, sample_gwas_df):
        """plot_manhattan_qq should work with bokeh backend."""
        pytest.importorskip("bokeh")

        plotter = LocusZoomPlotter(species="human", backend="bokeh")
        fig = plotter.plot_manhattan_qq(sample_gwas_df)
        assert fig is not None

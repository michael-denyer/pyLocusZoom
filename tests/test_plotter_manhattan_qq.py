"""Tests for Manhattan and QQ plot methods in ManhattanPlotter."""

import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import PathCollection, PolyCollection

from pylocuszoom import GenomeWideConfig
from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.manhattan_plotter import ManhattanPlotter
from tests.conftest import FIGURE_TYPES
from tests.figure_probes import PROBES

GENOMEWIDE_LINE = pytest.approx([-np.log10(5e-8)])


def _xtick_labels(ax):
    return [tick.get_text() for tick in ax.get_xticklabels()]


def _point_count(ax):
    return sum(
        len(c.get_offsets()) for c in ax.collections if isinstance(c, PathCollection)
    )


def _texts(ax):
    return [text.get_text() for text in ax.texts]


class TestPlotManhattan:
    """Tests for plot_manhattan method."""

    def test_plot_manhattan_with_custom_columns(self, manhattan_plotter):
        """Configured column names are the ones read: every row is drawn."""
        df = pd.DataFrame(
            {
                "chromosome": [1, 1, 2],
                "position": [1e6, 2e6, 1e6],
                "pvalue": [1e-8, 0.01, 0.5],
            }
        )
        fig = manhattan_plotter.plot_manhattan(
            df,
            config=GenomeWideConfig(
                chrom_col="chromosome", pos_col="position", p_col="pvalue"
            ),
        )
        assert _point_count(fig.get_axes()[0]) == 3

    def test_plot_manhattan_with_species_order(self):
        """Chromosomes are laid out in the species order, not the input order."""
        df = pd.DataFrame(
            {"chr": ["X", "2", "1"], "pos": [1e6, 1e6, 1e6], "p_value": [0.1, 0.2, 0.3]}
        )
        fig = ManhattanPlotter(species="canine").plot_manhattan(df)
        assert _xtick_labels(fig.get_axes()[0]) == ["1", "2", "X"]

    def test_plot_manhattan_with_custom_order(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        """A custom chromosome order replaces the species order."""
        fig = manhattan_plotter.plot_manhattan(
            manhattan_gwas_df,
            config=GenomeWideConfig(custom_chrom_order=["3", "2", "1"]),
        )
        assert _xtick_labels(fig.get_axes()[0]) == ["3", "2", "1"]

    def test_plot_manhattan_shows_significance_line(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        """plot_manhattan draws the genome-wide line at 5e-8 by default."""
        fig = manhattan_plotter.plot_manhattan(manhattan_gwas_df)
        assert PROBES["matplotlib"].hline_levels(fig) == GENOMEWIDE_LINE

    def test_plot_manhattan_custom_threshold(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        """A custom significance threshold moves the line."""
        fig = manhattan_plotter.plot_manhattan(
            manhattan_gwas_df, significance_threshold=1e-5
        )
        assert PROBES["matplotlib"].hline_levels(fig) == pytest.approx([5.0])

    def test_plot_manhattan_no_threshold(self, manhattan_plotter, manhattan_gwas_df):
        """significance_threshold=None draws no line."""
        fig = manhattan_plotter.plot_manhattan(
            manhattan_gwas_df, significance_threshold=None
        )
        assert PROBES["matplotlib"].hline_levels(fig) == []

    def test_plot_manhattan_with_figsize(self, manhattan_plotter, manhattan_gwas_df):
        """plot_manhattan should accept figsize parameter."""
        fig = manhattan_plotter.plot_manhattan(manhattan_gwas_df, figsize=(12, 4))
        assert fig.get_size_inches()[0] == pytest.approx(12, rel=0.1)
        assert fig.get_size_inches()[1] == pytest.approx(4, rel=0.1)

    def test_plot_manhattan_with_title(self, manhattan_plotter, manhattan_gwas_df):
        """plot_manhattan should accept title parameter."""
        fig = manhattan_plotter.plot_manhattan(
            manhattan_gwas_df, title="Test Manhattan"
        )
        ax = fig.get_axes()[0]
        assert "Test Manhattan" in ax.get_title()

    def test_plot_manhattan_validates_columns(self, manhattan_plotter):
        """plot_manhattan should raise on missing columns."""
        df = pd.DataFrame({"wrong": [1], "columns": [2]})
        with pytest.raises(ValueError, match="Missing columns"):
            manhattan_plotter.plot_manhattan(df)

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_plot_manhattan_on_every_backend(self, backend, manhattan_gwas_df):
        """Each backend draws one panel with the genome-wide line."""
        plotter = ManhattanPlotter(species="human", backend=backend)

        fig = plotter.plot_manhattan(manhattan_gwas_df)

        probe = PROBES[backend]
        assert isinstance(fig, FIGURE_TYPES[backend])
        assert probe.panel_count(fig) == 1
        assert probe.hline_levels(fig) == GENOMEWIDE_LINE


class TestPlotQQ:
    """Tests for plot_qq method."""

    @pytest.fixture
    def sample_pvalues_df(self):
        """Sample DataFrame with p-values for QQ plot."""
        rng = np.random.default_rng(42)
        return pd.DataFrame({"p_value": rng.uniform(0, 1, 1000)})

    @pytest.fixture
    def default_manhattan_plotter(self):
        """Create a plotter instance."""
        return ManhattanPlotter()

    @staticmethod
    def _bands(ax):
        return [
            c
            for c in ax.collections
            if isinstance(c, PolyCollection) and not isinstance(c, PathCollection)
        ]

    def test_plot_qq_with_custom_column(self, default_manhattan_plotter):
        """A configured p-value column is the one plotted: every value is drawn."""
        df = pd.DataFrame({"pvalue": np.random.default_rng(0).uniform(0, 1, 100)})
        fig = default_manhattan_plotter.plot_qq(
            df, config=GenomeWideConfig(p_col="pvalue")
        )
        assert _point_count(fig.get_axes()[0]) == 100

    def test_plot_qq_shows_confidence_band(
        self, default_manhattan_plotter, sample_pvalues_df
    ):
        """plot_qq shades one confidence band by default."""
        fig = default_manhattan_plotter.plot_qq(sample_pvalues_df)
        assert len(self._bands(fig.get_axes()[0])) == 1

    def test_plot_qq_no_confidence_band(
        self, default_manhattan_plotter, sample_pvalues_df
    ):
        """show_confidence_band=False shades nothing."""
        fig = default_manhattan_plotter.plot_qq(
            sample_pvalues_df, show_confidence_band=False
        )
        assert self._bands(fig.get_axes()[0]) == []

    def test_plot_qq_shows_diagonal(self, default_manhattan_plotter, sample_pvalues_df):
        """plot_qq draws the y = x reference line."""
        fig = default_manhattan_plotter.plot_qq(sample_pvalues_df)
        (line,) = fig.get_axes()[0].get_lines()
        assert list(line.get_xdata()) == list(line.get_ydata())

    def test_plot_qq_shows_lambda(self, default_manhattan_plotter, sample_pvalues_df):
        """show_lambda puts the genomic inflation factor in the title."""
        fig = default_manhattan_plotter.plot_qq(sample_pvalues_df, show_lambda=True)
        assert fig.get_axes()[0].get_title().startswith("QQ Plot (λ = ")

    def test_plot_qq_with_figsize(self, default_manhattan_plotter, sample_pvalues_df):
        """plot_qq should accept figsize parameter."""
        fig = default_manhattan_plotter.plot_qq(sample_pvalues_df, figsize=(6, 6))
        assert fig.get_size_inches()[0] == pytest.approx(6, rel=0.1)
        assert fig.get_size_inches()[1] == pytest.approx(6, rel=0.1)

    def test_plot_qq_with_title(self, default_manhattan_plotter, sample_pvalues_df):
        """plot_qq should accept title parameter."""
        fig = default_manhattan_plotter.plot_qq(sample_pvalues_df, title="Test QQ Plot")
        ax = fig.get_axes()[0]
        assert "Test QQ" in ax.get_title()

    def test_plot_qq_validates_columns(self, default_manhattan_plotter):
        """plot_qq should raise on missing p-value column."""
        df = pd.DataFrame({"wrong": [1, 2, 3]})
        with pytest.raises(ValueError, match="not found"):
            default_manhattan_plotter.plot_qq(df)

    def test_plot_qq_handles_all_nan(self, default_manhattan_plotter):
        """plot_qq should raise on all NaN p-values."""
        df = pd.DataFrame({"p_value": [np.nan, np.nan, np.nan]})
        with pytest.raises(ValueError, match="No valid"):
            default_manhattan_plotter.plot_qq(df)

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_plot_qq_on_every_backend(self, backend, sample_pvalues_df):
        """Each backend draws the QQ plot as one panel."""
        plotter = ManhattanPlotter(backend=backend)

        fig = plotter.plot_qq(sample_pvalues_df)

        assert isinstance(fig, FIGURE_TYPES[backend])
        assert PROBES[backend].panel_count(fig) == 1


class TestPlotManhattanCategorical:
    """Tests for categorical Manhattan plots (PheWAS-style)."""

    @pytest.fixture
    def sample_phewas_df(self):
        """Sample PheWAS-style DataFrame."""
        return pd.DataFrame(
            {
                "category": ["cardio", "cardio", "neuro", "neuro", "immuno"],
                "phenotype": ["BP", "HR", "AD", "PD", "RA"],
                "p_value": [1e-10, 0.01, 1e-6, 0.5, 1e-4],
            }
        )

    @pytest.fixture
    def default_manhattan_plotter(self):
        """Create a plotter instance."""
        return ManhattanPlotter()

    def test_plot_manhattan_categorical(
        self, default_manhattan_plotter, sample_phewas_df
    ):
        """category_col ticks the x-axis by category, alphabetically by default."""
        fig = default_manhattan_plotter.plot_manhattan(
            sample_phewas_df,
            category_col="category",
        )
        ax = fig.get_axes()[0]
        assert _xtick_labels(ax) == ["cardio", "immuno", "neuro"]
        assert _point_count(ax) == 5

    def test_plot_manhattan_categorical_custom_order(
        self, default_manhattan_plotter, sample_phewas_df
    ):
        """category_order sets the left-to-right category order."""
        fig = default_manhattan_plotter.plot_manhattan(
            sample_phewas_df,
            category_col="category",
            category_order=["neuro", "cardio", "immuno"],
        )
        assert _xtick_labels(fig.get_axes()[0]) == ["neuro", "cardio", "immuno"]


class TestPlotManhattanStacked:
    """Tests for stacked Manhattan plots."""

    @pytest.fixture
    def manhattan_gwas_dfs(self):
        """Three multi-chromosome GWAS frames in the chrom/pos/p schema."""
        rng = np.random.default_rng(42)
        dfs = []
        for i in range(3):
            n_variants = 50
            dfs.append(
                pd.DataFrame(
                    {
                        "chr": np.repeat([1, 2], [25, 25]),
                        "pos": np.concatenate(
                            [
                                np.sort(rng.integers(int(1e6), int(1e8), 25)),
                                np.sort(rng.integers(int(1e6), int(1e8), 25)),
                            ]
                        ),
                        "p_value": rng.uniform(1e-10, 1, n_variants),
                    }
                )
            )
        return dfs

    def test_plot_manhattan_stacked_creates_multiple_panels(
        self, manhattan_plotter, manhattan_gwas_dfs
    ):
        """plot_manhattan_stacked should create one panel per DataFrame."""
        fig = manhattan_plotter.plot_manhattan_stacked(manhattan_gwas_dfs)
        axes = fig.get_axes()
        # Should have 3 panels
        assert len(axes) == 3

    def test_plot_manhattan_stacked_with_panel_labels(
        self, manhattan_plotter, manhattan_gwas_dfs
    ):
        """Each panel carries its own label."""
        labels = ["Study A", "Study B", "Study C"]
        fig = manhattan_plotter.plot_manhattan_stacked(
            manhattan_gwas_dfs, panel_labels=labels
        )
        assert [_texts(ax) for ax in fig.get_axes()] == [[label] for label in labels]

    def test_plot_manhattan_stacked_validates_label_count(
        self, manhattan_plotter, manhattan_gwas_dfs
    ):
        """plot_manhattan_stacked should raise if panel_labels length mismatch."""
        with pytest.raises(ValueError, match="length"):
            manhattan_plotter.plot_manhattan_stacked(
                manhattan_gwas_dfs, panel_labels=["A", "B"]
            )

    def test_plot_manhattan_stacked_with_figsize(
        self, manhattan_plotter, manhattan_gwas_dfs
    ):
        """plot_manhattan_stacked should accept figsize parameter."""
        fig = manhattan_plotter.plot_manhattan_stacked(
            manhattan_gwas_dfs, figsize=(14, 10)
        )
        assert fig.get_size_inches()[0] == pytest.approx(14, rel=0.1)

    def test_plot_manhattan_stacked_single_df(
        self, manhattan_plotter, manhattan_gwas_dfs
    ):
        """A one-frame stack is a single panel."""
        fig = manhattan_plotter.plot_manhattan_stacked([manhattan_gwas_dfs[0]])
        assert len(fig.get_axes()) == 1

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_plot_manhattan_stacked_on_every_backend(self, backend, manhattan_gwas_dfs):
        """Each backend draws one panel per frame."""
        plotter = ManhattanPlotter(species="human", backend=backend)

        fig = plotter.plot_manhattan_stacked(manhattan_gwas_dfs)

        assert isinstance(fig, FIGURE_TYPES[backend])
        assert PROBES[backend].panel_count(fig) == 3


class TestPlotManhattanQQSideBySide:
    """Tests for side-by-side Manhattan and QQ plots."""

    def test_plot_manhattan_qq_creates_two_panels(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        """plot_manhattan_qq should create two side-by-side panels."""
        fig = manhattan_plotter.plot_manhattan_qq(manhattan_gwas_df)
        axes = fig.get_axes()
        # Should have 2 panels (Manhattan + QQ)
        assert len(axes) == 2

    def test_plot_manhattan_qq_with_title(self, manhattan_plotter, manhattan_gwas_df):
        """title becomes the figure title above both panels."""
        fig = manhattan_plotter.plot_manhattan_qq(
            manhattan_gwas_df, title="Combined Plot"
        )
        assert fig.get_suptitle() == "Combined Plot"

    def test_plot_manhattan_qq_with_figsize(self, manhattan_plotter, manhattan_gwas_df):
        """plot_manhattan_qq should accept figsize parameter."""
        fig = manhattan_plotter.plot_manhattan_qq(manhattan_gwas_df, figsize=(16, 5))
        assert fig.get_size_inches()[0] == pytest.approx(16, rel=0.1)

    def test_plot_manhattan_qq_custom_columns(self, manhattan_plotter):
        """Configured column names feed both panels: every row is drawn twice."""
        df = pd.DataFrame(
            {
                "chromosome": [1, 1, 2],
                "position": [1e6, 2e6, 1e6],
                "pvalue": [1e-8, 0.01, 0.5],
            }
        )
        fig = manhattan_plotter.plot_manhattan_qq(
            df,
            config=GenomeWideConfig(
                chrom_col="chromosome", pos_col="position", p_col="pvalue"
            ),
        )
        assert [_point_count(ax) for ax in fig.get_axes()] == [3, 3]

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_plot_manhattan_qq_on_every_backend(self, backend, manhattan_gwas_df):
        """Each backend draws the Manhattan and QQ panels."""
        plotter = ManhattanPlotter(species="human", backend=backend)

        fig = plotter.plot_manhattan_qq(manhattan_gwas_df)

        assert isinstance(fig, FIGURE_TYPES[backend])
        assert PROBES[backend].panel_count(fig) == 2


class TestPlotManhattanQQStacked:
    """Tests for plot_manhattan_qq_stacked method."""

    @pytest.fixture
    def manhattan_str_chrom_gwas_dfs(self):
        """Two GWAS frames in the chrom/pos/p schema with string chromosomes."""
        rng = np.random.default_rng(42)
        dfs = []
        for _ in range(2):
            data = []
            for chrom in [1, 2, 3]:
                n = 50
                positions = np.sort(rng.integers(int(1e6), int(5e7), n))
                pvalues = rng.uniform(0, 1, n)
                # Add some significant hits
                pvalues[:3] = [1e-10, 1e-8, 1e-6]
                for i in range(n):
                    data.append(
                        {"chr": str(chrom), "pos": positions[i], "p_value": pvalues[i]}
                    )
            dfs.append(pd.DataFrame(data))
        return dfs

    def test_plot_manhattan_qq_stacked_creates_correct_panels(
        self, manhattan_plotter, manhattan_str_chrom_gwas_dfs
    ):
        """plot_manhattan_qq_stacked should create n_gwas * 2 panels (Manhattan + QQ each)."""
        fig = manhattan_plotter.plot_manhattan_qq_stacked(manhattan_str_chrom_gwas_dfs)
        axes = fig.get_axes()
        # Should have 4 panels (2 GWAS * 2 plots each)
        assert len(axes) == 4

    def test_plot_manhattan_qq_stacked_with_panel_labels(
        self, manhattan_plotter, manhattan_str_chrom_gwas_dfs
    ):
        """Each study's Manhattan panel carries its label; QQ panels carry none."""
        fig = manhattan_plotter.plot_manhattan_qq_stacked(
            manhattan_str_chrom_gwas_dfs, panel_labels=["Study A", "Study B"]
        )
        assert [_texts(ax) for ax in fig.get_axes()] == [
            ["Study A"],
            [],
            ["Study B"],
            [],
        ]

    def test_plot_manhattan_qq_stacked_with_title(
        self, manhattan_plotter, manhattan_str_chrom_gwas_dfs
    ):
        """title becomes the figure title above the grid."""
        fig = manhattan_plotter.plot_manhattan_qq_stacked(
            manhattan_str_chrom_gwas_dfs, title="Multi-study GWAS"
        )
        assert fig.get_suptitle() == "Multi-study GWAS"

    def test_plot_manhattan_qq_stacked_with_figsize(
        self, manhattan_plotter, manhattan_str_chrom_gwas_dfs
    ):
        """plot_manhattan_qq_stacked should accept figsize parameter."""
        fig = manhattan_plotter.plot_manhattan_qq_stacked(
            manhattan_str_chrom_gwas_dfs, figsize=(16, 10)
        )
        assert fig.get_size_inches()[0] == pytest.approx(16, rel=0.1)

    def test_plot_manhattan_qq_stacked_three_studies(self, manhattan_plotter):
        """plot_manhattan_qq_stacked should work with three GWAS datasets."""
        dfs = []
        for _ in range(3):
            data = [
                {"chr": "1", "pos": 1e6, "p_value": 1e-8},
                {"chr": "1", "pos": 2e6, "p_value": 0.01},
                {"chr": "2", "pos": 1e6, "p_value": 0.5},
            ]
            dfs.append(pd.DataFrame(data))
        fig = manhattan_plotter.plot_manhattan_qq_stacked(dfs)
        axes = fig.get_axes()
        assert len(axes) == 6  # 3 GWAS * 2 plots each

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_plot_manhattan_qq_stacked_on_every_backend(
        self, backend, manhattan_str_chrom_gwas_dfs
    ):
        """Each backend draws a Manhattan and a QQ panel per study."""
        plotter = ManhattanPlotter(species="human", backend=backend)

        fig = plotter.plot_manhattan_qq_stacked(manhattan_str_chrom_gwas_dfs)

        assert isinstance(fig, FIGURE_TYPES[backend])
        assert PROBES[backend].panel_count(fig) == 4


class TestYlimClamp:
    """Regression: Manhattan variants must never emit a degenerate ylim(0, 0).

    When every p-value rounds to 1, ``-log10(p) == 0`` everywhere so ``y_max``
    would be 0 and the historical ``ylim(0, 0.0)`` collapsed the axis. The
    plotter now floors the upper bound at 1.0 via ``_padded_ymax``. Covers the
    four Manhattan entry points. (All-NaN p-values are rejected earlier by the
    p-value validator, so no ``_padded_ymax`` codepath sees ``NaN`` in
    practice.)
    """

    @staticmethod
    def _flat_df():
        return pd.DataFrame(
            {
                "chr": np.repeat([1, 2, 3], 10),
                "pos": np.tile(np.arange(int(1e6), int(1e6) + 10) * 1000, 3),
                "p_value": [1.0] * 30,
            }
        )

    def test_plot_manhattan_ylim_floor(self, manhattan_plotter):
        fig = manhattan_plotter.plot_manhattan(self._flat_df())
        assert fig.get_axes()[0].get_ylim()[1] >= 1.0

    def test_plot_manhattan_stacked_ylim_floor(self, manhattan_plotter):
        fig = manhattan_plotter.plot_manhattan_stacked(
            [self._flat_df(), self._flat_df()]
        )
        for ax in fig.get_axes():
            assert ax.get_ylim()[1] >= 1.0

    def test_plot_manhattan_qq_ylim_floor(self, manhattan_plotter):
        df = self._flat_df()
        # Keep one real p so the QQ half doesn't trip its own all-NaN guard.
        df.loc[0, "p_value"] = 0.5
        fig = manhattan_plotter.plot_manhattan_qq(df)
        # Manhattan axis is the wider one; QQ is roughly square.
        manhattan_ax = max(fig.get_axes(), key=lambda a: a.get_position().width)
        assert manhattan_ax.get_ylim()[1] >= 1.0

    def test_plot_manhattan_qq_stacked_ylim_floor(self, manhattan_plotter):
        dfs = []
        for _ in range(2):
            df = self._flat_df()
            df.loc[0, "p_value"] = 0.5
            dfs.append(df)
        fig = manhattan_plotter.plot_manhattan_qq_stacked(dfs)
        # Pick the widest axis per panel row; those are Manhattan panels.
        axes_by_y = sorted(fig.get_axes(), key=lambda a: -a.get_position().y0)
        for panel in range(2):
            row = axes_by_y[panel * 2 : panel * 2 + 2]
            manhattan_ax = max(row, key=lambda a: a.get_position().width)
            assert manhattan_ax.get_ylim()[1] >= 1.0


class TestYlimCoversThresholdLines:
    """The y-limit reaches the threshold lines when every point sits below them."""

    THRESHOLD = 5e-8

    @staticmethod
    def _weak_df():
        return pd.DataFrame(
            {
                "chr": np.repeat([1, 2, 3], 10),
                "pos": np.tile(np.arange(1, 11) * 1_000_000, 3),
                "p_value": np.linspace(1e-5, 1.0, 30),
            }
        )

    def test_plot_manhattan_shows_significance_line(self, manhattan_plotter):
        fig = manhattan_plotter.plot_manhattan(
            self._weak_df(), significance_threshold=self.THRESHOLD
        )

        assert fig.get_axes()[0].get_ylim()[1] > -np.log10(self.THRESHOLD)

    def test_plot_manhattan_qq_shows_suggestive_only_line(self, manhattan_plotter):
        fig = manhattan_plotter.plot_manhattan_qq(
            self._weak_df(), significance_threshold=None, suggestive_threshold=1e-6
        )

        assert fig.get_axes()[0].get_ylim()[1] > 6

    def test_plot_manhattan_stacked_shows_significance_line(self, manhattan_plotter):
        fig = manhattan_plotter.plot_manhattan_stacked(
            [self._weak_df(), self._weak_df()], significance_threshold=self.THRESHOLD
        )

        for ax in fig.get_axes():
            assert ax.get_ylim()[1] > -np.log10(self.THRESHOLD)

    def test_plot_manhattan_categorical_shows_significance_line(
        self, manhattan_plotter
    ):
        df = pd.DataFrame({"cat": ["a", "b", "c"], "p_value": [0.1, 0.01, 1e-4]})

        fig = manhattan_plotter.plot_manhattan(
            df, category_col="cat", significance_threshold=self.THRESHOLD
        )

        assert fig.get_axes()[0].get_ylim()[1] > -np.log10(self.THRESHOLD)

    def test_miami_bottom_panel_stays_inverted(self):
        from pylocuszoom import MiamiPlotter

        fig = MiamiPlotter(species="human").plot_miami(
            self._weak_df(),
            self._weak_df(),
            top_threshold=self.THRESHOLD,
            bottom_threshold=self.THRESHOLD,
        )

        top, bottom = fig.get_axes()
        line = -np.log10(self.THRESHOLD)
        assert top.get_ylim()[1] > line
        assert bottom.get_ylim()[0] > line
        assert bottom.get_ylim()[1] == 0


class TestPlotManhattanQQOptions:
    """Suggestive line, caller-supplied lambda and footer on plot_manhattan_qq."""

    @staticmethod
    def _hline_levels(ax):
        return sorted(round(line.get_ydata()[0], 6) for line in ax.get_lines())

    def test_defaults_draw_one_line_and_computed_lambda(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        fig = manhattan_plotter.plot_manhattan_qq(manhattan_gwas_df)

        manhattan_ax, qq_ax = fig.get_axes()
        assert self._hline_levels(manhattan_ax) == [round(-np.log10(5e-8), 6)]
        assert qq_ax.get_title().startswith("QQ Plot (λ = ")
        assert fig.texts == []

    def test_suggestive_line_is_drawn_below_the_genomewide_line(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        fig = manhattan_plotter.plot_manhattan_qq(
            manhattan_gwas_df, suggestive_threshold=1e-5
        )

        assert self._hline_levels(fig.get_axes()[0]) == [5.0, round(-np.log10(5e-8), 6)]

    def test_caller_lambda_replaces_the_computed_one(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        fig = manhattan_plotter.plot_manhattan_qq(manhattan_gwas_df, lambda_gc=1.2345)

        assert fig.get_axes()[1].get_title() == "QQ Plot (λ = 1.234)"

    def test_footer_is_written_under_the_panels(
        self, manhattan_plotter, manhattan_gwas_df
    ):
        fig = manhattan_plotter.plot_manhattan_qq(
            manhattan_gwas_df, footer="n = 1,234 dogs"
        )

        assert [text.get_text() for text in fig.texts] == ["n = 1,234 dogs"]

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_options_render_on_every_backend(self, backend, manhattan_gwas_df):
        plotter = ManhattanPlotter(species="human", backend=backend)

        fig = plotter.plot_manhattan_qq(
            manhattan_gwas_df,
            suggestive_threshold=1e-5,
            lambda_gc=1.05,
            footer="footer text",
        )

        assert isinstance(fig, FIGURE_TYPES[backend])
        if backend == "plotly":
            texts = [a.text for a in fig.layout.annotations]
            assert "<i>footer text</i>" in texts
        elif backend == "bokeh":
            from bokeh.models import Div

            assert any(
                isinstance(model, Div) and model.text == "footer text"
                for model in fig.select({"type": Div})
            )

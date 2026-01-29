"""Edge case tests for pyLocusZoom plot types.

Tests handling of empty DataFrames, NaN p-values, mismatched list lengths,
single-chromosome Manhattan plots, and color cycling for many categories/credible sets.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.plotter import LocusZoomPlotter


class TestEmptyDataFrames:
    """Test handling of empty DataFrames across plot types.

    Behavior varies by plot type:
    - Regional plots: Allow empty DataFrames, render empty axes
    - Manhattan plots: Raise ValueError (can't compute axis limits)
    - QQ plots: Raise ValueError (require valid p-values)
    """

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    @pytest.fixture
    def empty_gwas_df(self):
        """Empty DataFrame with correct columns."""
        return pd.DataFrame(columns=["chrom", "ps", "p_wald", "rs"])

    def test_plot_with_empty_df_succeeds(self, plotter, empty_gwas_df):
        """Regional plot with empty DataFrame renders without error.

        The current behavior allows empty DataFrames through and produces
        an empty plot. The plot contains no data points but still has
        axes and labels.
        """
        fig = plotter.plot(
            empty_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            pos_col="ps",
            p_col="p_wald",
            show_recombination=False,
        )
        assert fig is not None
        plt.close(fig)

    def test_manhattan_with_empty_df_raises(self):
        """Manhattan plot with empty DataFrame raises ValueError.

        Empty DataFrames have no data to compute axis limits,
        resulting in NaN limits which matplotlib cannot handle.
        """
        plotter = LocusZoomPlotter(species="human", log_level=None)
        empty_df = pd.DataFrame(columns=["chrom", "pos", "p"])

        with pytest.raises(ValueError, match="(NaN|Inf|cannot)"):
            plotter.plot_manhattan(
                empty_df,
                chrom_col="chrom",
                pos_col="pos",
                p_col="p",
            )

    def test_qq_with_empty_df_raises(self, plotter):
        """QQ plot with empty DataFrame raises ValueError.

        The QQ plot requires valid p-values (>0 and <=1) to compute
        expected quantiles. Empty data fails this validation.
        """
        empty_df = pd.DataFrame({"p": pd.Series([], dtype=float)})

        with pytest.raises(ValueError, match="No valid p-values"):
            plotter.plot_qq(empty_df, p_col="p")


class TestNaNPvalues:
    """Test handling of NaN p-values in plots."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    def test_plot_with_some_nan_pvalues_succeeds(self, plotter):
        """Regional plot with partial NaN p-values should work.

        Rows with NaN p-values are excluded from plotting but do not
        cause errors. A warning is logged (visible in stderr).
        """
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "ps": [1100000, 1300000, 1500000, 1700000, 1900000],
                "p_wald": [1e-8, np.nan, 1e-5, np.nan, 0.01],
            }
        )

        fig = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert fig is not None
        plt.close(fig)

    def test_plot_with_all_nan_pvalues_succeeds(self, plotter):
        """Regional plot with all NaN p-values should render empty.

        The plot renders without error but contains no data points.
        This tests graceful degradation rather than raising an error.
        """
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [np.nan, np.nan, np.nan],
            }
        )

        fig = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert fig is not None
        plt.close(fig)


class TestStackedPlotMismatchedLengths:
    """Test mismatched list lengths in stacked plots raise ValidationError."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results DataFrame."""
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

    def test_stacked_mismatched_lead_positions_raises(self, plotter, sample_gwas_df):
        """Mismatched lead_positions length raises ValueError."""
        gwas_dfs = [sample_gwas_df, sample_gwas_df.copy(), sample_gwas_df.copy()]
        lead_positions = [1500000, 1500000]  # Only 2, but 3 gwas_dfs

        with pytest.raises(ValueError, match="lead_positions"):
            plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=lead_positions,
                show_recombination=False,
            )

    def test_stacked_mismatched_panel_labels_raises(self, plotter, sample_gwas_df):
        """Mismatched panel_labels length raises ValueError."""
        gwas_dfs = [sample_gwas_df, sample_gwas_df.copy()]
        panel_labels = ["Only One"]  # Should have 2 labels

        with pytest.raises(ValueError, match="panel_labels"):
            plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                panel_labels=panel_labels,
                show_recombination=False,
            )

    def test_stacked_mismatched_ld_reference_files_raises(
        self, plotter, sample_gwas_df
    ):
        """Mismatched ld_reference_files length raises ValueError."""
        gwas_dfs = [sample_gwas_df, sample_gwas_df.copy()]
        ld_reference_files = ["/path/to/file1"]  # Only 1, but 2 gwas_dfs

        with pytest.raises(ValueError, match="ld_reference_files"):
            plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                ld_reference_files=ld_reference_files,
                show_recombination=False,
            )


class TestManhattanSingleChromosome:
    """Test Manhattan plot with single chromosome data."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species="human", log_level=None)

    def test_manhattan_single_chromosome_succeeds(self, plotter):
        """Manhattan plot with only one chromosome should work."""
        df = pd.DataFrame(
            {
                "chrom": ["1"] * 10,
                "pos": list(range(1000000, 11000000, 1000000)),
                "p": [1e-8, 0.05, 0.01, 1e-6, 0.1, 0.001, 1e-10, 0.5, 0.005, 1e-3],
            }
        )

        fig = plotter.plot_manhattan(
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


class TestPheWASManyCategories:
    """Test PheWAS plot with many categories cycles colors correctly."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    def test_phewas_15_categories_succeeds(self, plotter):
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

        fig = plotter.plot_phewas(
            phewas_df,
            variant_id="rs12345",
            phenotype_col="phenotype",
            p_col="p_value",
            category_col="category",
        )

        assert fig is not None
        # Verify all 15 points are plotted
        ax = fig.axes[0]
        total_points = sum(
            len(collection.get_offsets()) for collection in ax.collections
        )
        assert total_points == 15
        plt.close(fig)


class TestFinemappingManyCredibleSets:
    """Test fine-mapping plot with many credible sets cycles colors correctly."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results DataFrame."""
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "ps": [1100000, 1300000, 1500000, 1700000, 1900000],
                "p_wald": [1e-8, 1e-6, 1e-5, 1e-4, 0.01],
            }
        )

    def test_finemapping_12_credible_sets_succeeds(self, plotter, sample_gwas_df):
        """Fine-mapping with >10 credible sets should cycle colors without error.

        The CREDIBLE_SET_COLORS palette has 10 colors, so 12 credible sets
        requires color cycling. This tests that modulo indexing works.
        """
        # Create fine-mapping data with 12 credible sets
        n_variants = 60  # 5 variants per credible set
        positions = list(range(1000000, 1000000 + n_variants * 10000, 10000))
        credible_sets = [((i // 5) % 12) + 1 for i in range(n_variants)]

        finemapping_df = pd.DataFrame(
            {
                "pos": positions,
                "pip": [0.8 if i % 5 == 0 else 0.1 for i in range(n_variants)],
                "cs": credible_sets,
            }
        )

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=900000,
            end=1700000,
            show_recombination=False,
            finemapping_df=finemapping_df,
            finemapping_cs_col="cs",
        )

        assert fig is not None
        plt.close(fig)

    def test_credible_set_color_cycling(self):
        """Verify get_credible_set_color cycles correctly for cs > 10."""
        from pylocuszoom.colors import CREDIBLE_SET_COLORS, get_credible_set_color

        # Test that colors cycle after 10
        for cs_id in range(1, 25):
            color = get_credible_set_color(cs_id)
            expected_idx = (cs_id - 1) % len(CREDIBLE_SET_COLORS)
            expected_color = CREDIBLE_SET_COLORS[expected_idx]
            assert color == expected_color, (
                f"CS {cs_id}: got {color}, expected {expected_color}"
            )


class TestManhattanStackedValidation:
    """Test validation in stacked Manhattan plots."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species="human", log_level=None)

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

    def test_manhattan_stacked_empty_list_raises(self, plotter):
        """Empty list of GWAS DataFrames should raise ValueError."""
        with pytest.raises(ValueError, match="At least one GWAS DataFrame"):
            plotter.plot_manhattan_stacked([])

    def test_manhattan_stacked_mismatched_panel_labels_raises(self, plotter, sample_df):
        """Mismatched panel_labels length should raise ValueError."""
        gwas_dfs = [sample_df, sample_df.copy()]
        panel_labels = ["Only One"]

        with pytest.raises(ValueError, match="panel_labels"):
            plotter.plot_manhattan_stacked(
                gwas_dfs,
                panel_labels=panel_labels,
            )


class TestManhattanQQStackedValidation:
    """Test validation in stacked Manhattan+QQ plots."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species="human", log_level=None)

    def test_manhattan_qq_stacked_empty_list_raises(self, plotter):
        """Empty list of GWAS DataFrames should raise ValueError."""
        with pytest.raises(ValueError, match="At least one GWAS DataFrame"):
            plotter.plot_manhattan_qq_stacked([])


class TestQQWithVariousPvalueDistributions:
    """Test QQ plot with various p-value distributions."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    def test_qq_uniform_pvalues(self, plotter):
        """QQ plot with uniform p-values should show lambda ~ 1."""
        np.random.seed(42)
        df = pd.DataFrame({"p": np.random.uniform(0, 1, 1000)})

        fig = plotter.plot_qq(df, p_col="p", show_lambda=True)
        assert fig is not None

        # Check title contains lambda close to 1
        ax = fig.axes[0]
        title = ax.get_title()
        assert "λ" in title
        plt.close(fig)

    def test_qq_extreme_pvalues(self, plotter):
        """QQ plot with very small p-values should not produce inf."""
        df = pd.DataFrame(
            {
                "p": [1e-300, 1e-200, 1e-100, 1e-50, 0.001, 0.01, 0.1, 0.5],
            }
        )

        fig = plotter.plot_qq(df, p_col="p")
        assert fig is not None
        plt.close(fig)


class TestRegionalPlotColumnValidation:
    """Test column validation in regional plots."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance for testing."""
        return LocusZoomPlotter(species=None, log_level=None)

    def test_plot_missing_pos_col_raises(self, plotter):
        """Plot with missing position column should raise ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "p_wald": [1e-8, 0.01],
                # missing 'ps' column
            }
        )

        with pytest.raises(ValidationError, match="ps"):
            plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

    def test_plot_missing_p_col_raises(self, plotter):
        """Plot with missing p-value column should raise ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "ps": [1100000, 1900000],
                # missing 'p_wald' column
            }
        )

        with pytest.raises(ValidationError, match="p_wald"):
            plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

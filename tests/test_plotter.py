"""Tests for LocusZoomPlotter class."""

import warnings
from pathlib import Path
from unittest.mock import MagicMock, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend
from pylocuszoom.backends.plotly_backend import PlotlyBackend
from pylocuszoom.plotter import LocusZoomPlotter


class TestBackendIntegration:
    """Tests for backend protocol integration."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results DataFrame."""
        np.random.seed(42)
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

    def test_plot_uses_backend_create_figure(self, sample_gwas_df):
        """plot() should use self._backend.create_figure() instead of plt.subplots()."""
        plotter = LocusZoomPlotter(species="canine")

        # Spy on the backend's create_figure method
        original_create_figure = plotter._backend.create_figure
        plotter._backend.create_figure = MagicMock(side_effect=original_create_figure)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )

        # Backend's create_figure should have been called
        plotter._backend.create_figure.assert_called()
        plt.close(fig)

    def test_plot_stacked_uses_backend_create_figure(self, sample_gwas_df):
        """plot_stacked() should use self._backend.create_figure()."""
        plotter = LocusZoomPlotter(species="canine")

        original_create_figure = plotter._backend.create_figure
        plotter._backend.create_figure = MagicMock(side_effect=original_create_figure)

        fig = plotter.plot_stacked(
            [sample_gwas_df, sample_gwas_df.copy()],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )

        plotter._backend.create_figure.assert_called()
        plt.close(fig)

    def test_default_backend_is_matplotlib(self):
        """Default backend should be matplotlib."""
        plotter = LocusZoomPlotter()
        assert isinstance(plotter._backend, MatplotlibBackend)

    def test_explicit_matplotlib_backend(self):
        """backend='matplotlib' should use MatplotlibBackend."""
        plotter = LocusZoomPlotter(backend="matplotlib")
        assert isinstance(plotter._backend, MatplotlibBackend)

    def test_explicit_plotly_backend(self):
        """backend='plotly' should use PlotlyBackend."""
        plotter = LocusZoomPlotter(backend="plotly")
        assert isinstance(plotter._backend, PlotlyBackend)

    def test_plotly_backend_creates_figure(self, sample_gwas_df):
        """plot() with backend='plotly' should create a plotly figure."""
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(species="canine", backend="plotly")

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )

        assert isinstance(fig, go.Figure)

    def test_plot_uses_backend_scatter(self, sample_gwas_df):
        """plot() should use self._backend.scatter() for association points."""
        plotter = LocusZoomPlotter(species="canine")

        original_scatter = plotter._backend.scatter
        plotter._backend.scatter = MagicMock(side_effect=original_scatter)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )

        plotter._backend.scatter.assert_called()
        plt.close(fig)

    def test_plot_uses_backend_axhline(self, sample_gwas_df):
        """plot() should use self._backend.axhline() for significance line."""
        plotter = LocusZoomPlotter(species="canine")

        original_axhline = plotter._backend.axhline
        plotter._backend.axhline = MagicMock(side_effect=original_axhline)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )

        plotter._backend.axhline.assert_called()
        plt.close(fig)

    def test_plot_uses_backend_axis_methods(self, sample_gwas_df):
        """plot() should use backend methods for axis configuration."""
        plotter = LocusZoomPlotter(species="canine")

        original_set_ylabel = plotter._backend.set_ylabel
        original_set_xlim = plotter._backend.set_xlim
        plotter._backend.set_ylabel = MagicMock(side_effect=original_set_ylabel)
        plotter._backend.set_xlim = MagicMock(side_effect=original_set_xlim)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )

        plotter._backend.set_ylabel.assert_called()
        plotter._backend.set_xlim.assert_called()
        plt.close(fig)


class TestLocusZoomPlotterInit:
    """Tests for LocusZoomPlotter initialization."""

    def test_default_species_is_canine(self):
        """Default species should be canine."""
        plotter = LocusZoomPlotter()
        assert plotter.species == "canine"

    def test_custom_species(self):
        """Should accept custom species."""
        plotter = LocusZoomPlotter(species="feline")
        assert plotter.species == "feline"

    def test_custom_plink_path(self):
        """Should accept custom PLINK path."""
        plotter = LocusZoomPlotter(plink_path="/custom/plink")
        assert plotter.plink_path == "/custom/plink"

    def test_custom_threshold(self):
        """Should accept custom genomewide threshold."""
        plotter = LocusZoomPlotter(genomewide_threshold=5e-8)
        assert plotter.genomewide_threshold == 5e-8

    def test_auto_genes_default_false(self):
        """auto_genes should be False by default for backward compatibility."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        assert plotter._auto_genes is False

    def test_auto_genes_can_be_enabled(self):
        """auto_genes=True should be accepted."""
        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)
        assert plotter._auto_genes is True


class TestAutoGenes:
    """Tests for automatic gene fetching from Ensembl."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS DataFrame for testing."""
        return pd.DataFrame(
            {
                "ps": [1100000, 1200000, 1300000, 1400000, 1500000],
                "p_wald": [1e-8, 1e-6, 1e-5, 1e-4, 0.01],
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            }
        )

    @pytest.fixture
    def sample_genes_df(self):
        """Sample gene DataFrame for testing."""
        return pd.DataFrame(
            {
                "chr": ["1", "1"],
                "start": [1100000, 1300000],
                "end": [1200000, 1400000],
                "gene_name": ["GENE1", "GENE2"],
                "strand": ["+", "-"],
            }
        )

    def test_plot_with_auto_genes_enabled(self, sample_gwas_df):
        """Test that auto_genes=True fetches genes from Ensembl."""
        # Mock the Ensembl API response
        mock_genes = pd.DataFrame(
            {
                "chr": ["1", "1"],
                "start": [1000000, 1500000],
                "end": [1200000, 1700000],
                "gene_name": ["GENE1", "GENE2"],
                "strand": ["+", "-"],
            }
        )

        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)

        with patch("pylocuszoom.plotter.get_genes_for_region", return_value=mock_genes):
            fig = plotter.plot(
                sample_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

        assert fig is not None

    def test_plot_auto_genes_disabled_by_default(self, sample_gwas_df):
        """Test that auto_genes=False by default (backward compatible)."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        # Should work without genes_df and without calling Ensembl
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
        )

        assert fig is not None

    def test_plot_auto_genes_respects_explicit_genes_df(
        self, sample_gwas_df, sample_genes_df
    ):
        """Test that explicit genes_df is used even when auto_genes=True."""
        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)

        with patch("pylocuszoom.plotter.get_genes_for_region") as mock_fetch:
            fig = plotter.plot(
                sample_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                genes_df=sample_genes_df,
            )

            # Ensembl should NOT be called when genes_df is provided
            mock_fetch.assert_not_called()

        assert fig is not None


class TestLocusZoomPlotterPlot:
    """Tests for LocusZoomPlotter.plot() method."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance."""
        return LocusZoomPlotter(species="canine")

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results DataFrame."""
        np.random.seed(42)
        n_snps = 50
        positions = np.sort(np.random.randint(1000000, 2000000, n_snps))
        return pd.DataFrame(
            {
                "rs": [f"rs{i}" for i in range(n_snps)],
                "chr": [1] * n_snps,
                "ps": positions,
                "p_wald": np.random.uniform(1e-10, 1, n_snps),
            }
        )

    @pytest.fixture
    def sample_genes_df(self):
        """Sample gene annotations."""
        return pd.DataFrame(
            {
                "chr": [1, 1, 1],
                "start": [1100000, 1400000, 1700000],
                "end": [1150000, 1500000, 1800000],
                "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
                "strand": ["+", "-", "+"],
            }
        )

    def test_creates_figure(self, plotter, sample_gwas_df):
        """Should create a matplotlib figure."""
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plots_with_gene_track(self, plotter, sample_gwas_df, sample_genes_df):
        """Should create plot with gene track when genes_df provided."""
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            genes_df=sample_genes_df,
        )
        # Should have 2 axes (association + gene track)
        assert len(fig.axes) >= 2
        plt.close(fig)

    def test_highlights_lead_snp(self, plotter, sample_gwas_df):
        """Should highlight lead SNP when lead_pos provided."""
        lead_pos = sample_gwas_df["ps"].iloc[0]
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            lead_pos=lead_pos,
        )
        plt.close(fig)
        # Test passes if no exception

    def test_handles_empty_dataframe(self, plotter):
        """Should handle empty GWAS DataFrame."""
        empty_df = pd.DataFrame(columns=["rs", "chr", "ps", "p_wald"])
        fig = plotter.plot(
            empty_df,
            chrom=1,
            start=1000000,
            end=2000000,
        )
        plt.close(fig)

    def test_custom_column_names(self, plotter):
        """Should work with custom column names."""
        df = pd.DataFrame(
            {
                "snp_id": ["rs1", "rs2", "rs3"],
                "position": [1100000, 1500000, 1900000],
                "pvalue": [1e-8, 1e-5, 1e-3],
            }
        )
        fig = plotter.plot(
            df,
            chrom=1,
            start=1000000,
            end=2000000,
            pos_col="position",
            p_col="pvalue",
            rs_col="snp_id",
        )
        plt.close(fig)

    def test_with_precomputed_ld(self, plotter, sample_gwas_df):
        """Should use pre-computed LD column when provided."""
        df = sample_gwas_df.copy()
        df["R2"] = np.random.uniform(0, 1, len(df))

        fig = plotter.plot(
            df,
            chrom=1,
            start=1000000,
            end=2000000,
            ld_col="R2",
        )
        plt.close(fig)

    def test_with_recombination_data(self, plotter, sample_gwas_df):
        """Should plot with recombination overlay when provided."""
        recomb_df = pd.DataFrame(
            {
                "pos": [1000000, 1200000, 1400000, 1600000, 1800000, 2000000],
                "rate": [0.5, 1.2, 2.5, 1.8, 0.8, 0.3],
            }
        )
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            recomb_df=recomb_df,
        )
        plt.close(fig)

    def test_disables_snp_labels(self, plotter, sample_gwas_df):
        """Should not add labels when snp_labels=False."""
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            snp_labels=False,
        )
        plt.close(fig)

    def test_disables_recombination(self, plotter, sample_gwas_df):
        """Should not show recombination when show_recombination=False."""
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        plt.close(fig)


class TestLocusZoomPlotterLdCalculation:
    """Tests for LD calculation integration."""

    @pytest.fixture
    def plotter(self):
        """Create plotter with mocked PLINK."""
        return LocusZoomPlotter(species="canine", plink_path="/mock/plink")

    def test_calculates_ld_when_reference_provided(self, plotter):
        """Should attempt LD calculation when ld_reference_file provided."""
        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom.plotter.calculate_ld") as mock_ld:
            mock_ld.return_value = pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "R2": [1.0, 0.8, 0.5],
                }
            )

            fig = plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_pos=1100000,
                ld_reference_file="/path/to/genotypes",
            )

            mock_ld.assert_called_once()
            plt.close(fig)


class TestLocusZoomPlotterRecombination:
    """Tests for recombination data handling."""

    def test_caches_recombination_data(self):
        """Should cache recombination data for repeated calls."""
        plotter = LocusZoomPlotter(species=None)  # No auto-download

        recomb_df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "rate": [0.5, 1.0, 0.5],
            }
        )

        # First call - no cache
        assert plotter._recomb_cache == {}

        # Manually add to cache (key includes genome_build)
        plotter._recomb_cache[(1, 1000000, 2000000, plotter.genome_build)] = recomb_df

        # Should return cached data
        result = plotter._get_recomb_for_region(1, 1000000, 2000000)
        assert result is not None
        assert len(result) == 3

    def test_recombination_overlay_does_not_distort_primary_ylim(self):
        """Primary y-axis limits should be unchanged when recombination is enabled.

        Regression test: recombination overlay was being plotted on the primary axis
        instead of a twin axis, causing GWAS y-limits to be rescaled by recomb rates.
        """
        plotter = LocusZoomPlotter(species=None)

        gwas_df = pd.DataFrame(
            {
                "rs": [f"rs{i}" for i in range(10)],
                "chr": [1] * 10,
                "ps": list(range(1000000, 2000000, 100000)),
                "p_wald": [1e-8, 1e-6, 1e-5, 1e-4, 0.01, 0.05, 0.1, 0.5, 0.8, 0.99],
            }
        )

        recomb_df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "rate": [50.0, 100.0, 75.0],  # High rates that would distort y-axis
            }
        )

        # Plot without recombination
        fig_no_recomb = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        ax_no_recomb = fig_no_recomb.axes[0]
        ylim_no_recomb = ax_no_recomb.get_ylim()

        # Plot with recombination
        fig_with_recomb = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            recomb_df=recomb_df,
        )
        ax_with_recomb = fig_with_recomb.axes[0]
        ylim_with_recomb = ax_with_recomb.get_ylim()

        # Primary y-axis limits should be the same
        assert ylim_no_recomb == ylim_with_recomb, (
            f"Recombination overlay distorted primary y-axis: "
            f"without={ylim_no_recomb}, with={ylim_with_recomb}"
        )

        plt.close(fig_no_recomb)
        plt.close(fig_with_recomb)


class TestPlotEdgeCases:
    """Tests for plot() edge cases and error handling."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance."""
        return LocusZoomPlotter(species="canine", plink_path="/mock/plink")

    def test_plot_raises_keyerror_when_rs_col_missing_with_ld_reference(self, plotter):
        """Bug: plot() should handle missing rs_col when ld_reference_file provided.

        Currently raises KeyError at line 264 when rs_col column doesn't exist
        but ld_reference_file is provided. Should either:
        1. Validate rs_col exists upfront and raise clear error, or
        2. Skip LD calculation gracefully with a warning
        """
        # GWAS data WITHOUT rs column
        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom.plotter.calculate_ld") as mock_ld:
            mock_ld.return_value = pd.DataFrame({"SNP": [], "R2": []})

            # This should NOT raise KeyError - should handle gracefully
            # Currently fails with: KeyError: 'rs'
            fig = plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_pos=1500000,
                ld_reference_file="/path/to/genotypes",
            )
            plt.close(fig)


class TestPlotStackedEdgeCases:
    """Tests for plot_stacked() edge cases and error handling."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance."""
        return LocusZoomPlotter(species="canine", log_level=None)

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

    def test_plot_stacked_validates_eqtl_columns(self, plotter, sample_gwas_df):
        """Bug: plot_stacked() should validate eQTL DataFrame has required columns.

        Currently bypasses validate_eqtl_df() and directly accesses 'pos' and
        'p_value' columns at lines 945-952, causing cryptic KeyError instead
        of helpful validation message.
        """
        from pylocuszoom.eqtl import EQTLValidationError

        # eQTL data with wrong column names
        bad_eqtl_df = pd.DataFrame(
            {
                "position": [1500000],  # Should be 'pos'
                "pval": [1e-6],  # Should be 'p_value'
            }
        )

        # Should raise EQTLValidationError with helpful message
        # Currently raises KeyError: 'pos'
        with pytest.raises(EQTLValidationError):
            plotter.plot_stacked(
                [sample_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=False,
                eqtl_df=bad_eqtl_df,
            )

    def test_plot_stacked_validates_list_lengths(self, plotter, sample_gwas_df):
        """Bug: plot_stacked() should error when list lengths don't match.

        Currently uses zip() which silently truncates the longer list.
        If user provides 3 GWAS DataFrames but only 2 lead_positions,
        the third GWAS is plotted without a lead SNP - confusing behavior.
        """
        gwas_dfs = [sample_gwas_df, sample_gwas_df.copy(), sample_gwas_df.copy()]
        lead_positions = [1500000, 1500000]  # Only 2, but 3 gwas_dfs

        # Should raise ValueError about mismatched lengths
        # Currently silently truncates - third GWAS has no lead SNP
        with pytest.raises(ValueError, match="lead_positions"):
            plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=lead_positions,
                show_recombination=False,
            )

    def test_plot_stacked_validates_panel_labels_length(self, plotter, sample_gwas_df):
        """Bug: panel_labels length should match gwas_dfs length."""
        gwas_dfs = [sample_gwas_df, sample_gwas_df.copy()]
        panel_labels = ["Only One"]  # Should have 2 labels

        # Should raise ValueError about mismatched lengths
        # Currently silently ignores - second panel has no label
        with pytest.raises(ValueError, match="panel_labels"):
            plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                panel_labels=panel_labels,
                show_recombination=False,
            )


class TestPValueValidation:
    """Tests for p-value validation and NaN handling."""

    @pytest.fixture
    def plotter(self):
        """Create plotter instance."""
        return LocusZoomPlotter(species=None)

    def test_plot_handles_nan_pvalues_with_warning(self, plotter):
        """Plot should handle NaN p-values and log a warning."""
        import numpy as np

        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, np.nan, 0.01],  # One NaN p-value
            }
        )

        # Should not raise, but should warn (captured by logging)
        fig = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        plt.close(fig)

    def test_plot_stacked_handles_all_nan_pvalues(self, plotter):
        """plot_stacked should handle region with all NaN p-values.

        Regression test: idxmin() on all-NaN series returns NaN,
        causing subsequent loc to fail.
        """
        import numpy as np

        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [np.nan, np.nan, np.nan],  # All NaN
            }
        )

        # Should not raise - should handle gracefully
        fig = plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        plt.close(fig)

    def test_plot_handles_out_of_range_pvalues(self, plotter):
        """Plot should handle p-values outside [0, 1] range."""
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [-0.1, 1.5, 0.05],  # Out of range values
            }
        )

        # Should not raise, but should warn
        fig = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        plt.close(fig)


class TestBackendEQTLFinemapping:
    """Tests for eQTL and fine-mapping support across all backends."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results DataFrame."""
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "ps": [1100000, 1300000, 1500000, 1700000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3, 0.01, 0.1],
            }
        )

    @pytest.fixture
    def sample_eqtl_df(self):
        """Sample eQTL DataFrame with effect sizes."""
        return pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],
                "p_value": [1e-6, 1e-4, 0.01],
                "gene": ["GENE1", "GENE1", "GENE1"],
                "effect_size": [0.5, -0.3, 0.8],
            }
        )

    @pytest.fixture
    def sample_eqtl_df_no_effect(self):
        """Sample eQTL DataFrame without effect sizes."""
        return pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],
                "p_value": [1e-6, 1e-4, 0.01],
                "gene": ["GENE1", "GENE1", "GENE1"],
            }
        )

    @pytest.fixture
    def sample_finemapping_df(self):
        """Sample fine-mapping DataFrame with credible sets."""
        return pd.DataFrame(
            {
                "pos": [1100000, 1300000, 1500000, 1700000, 1900000],
                "pip": [0.85, 0.12, 0.02, 0.45, 0.01],
                "cs": [1, 1, 0, 2, 0],
            }
        )

    def test_matplotlib_eqtl_with_effects(self, sample_gwas_df, sample_eqtl_df):
        """Matplotlib backend should handle eQTL panel with effect sizes."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE1",
        )

        assert fig is not None
        plt.close(fig)

    def test_matplotlib_eqtl_without_effects(
        self, sample_gwas_df, sample_eqtl_df_no_effect
    ):
        """Matplotlib backend should handle eQTL panel without effect sizes."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df_no_effect,
            eqtl_gene="GENE1",
        )

        assert fig is not None
        plt.close(fig)

    def test_matplotlib_finemapping(self, sample_gwas_df, sample_finemapping_df):
        """Matplotlib backend should handle fine-mapping panel."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
        )

        assert fig is not None
        plt.close(fig)

    def test_plotly_eqtl_with_effects(self, sample_gwas_df, sample_eqtl_df):
        """Plotly backend should handle eQTL panel with effect sizes without error."""
        plotter = LocusZoomPlotter(species=None, backend="plotly", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE1",
        )

        assert fig is not None
        # Plotly figures are go.Figure objects

    def test_plotly_eqtl_without_effects(
        self, sample_gwas_df, sample_eqtl_df_no_effect
    ):
        """Plotly backend should handle eQTL panel without effect sizes."""
        plotter = LocusZoomPlotter(species=None, backend="plotly", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df_no_effect,
            eqtl_gene="GENE1",
        )

        assert fig is not None

    def test_plotly_finemapping(self, sample_gwas_df, sample_finemapping_df):
        """Plotly backend should handle fine-mapping panel without error."""
        plotter = LocusZoomPlotter(species=None, backend="plotly", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
        )

        assert fig is not None

    def test_bokeh_eqtl_with_effects(self, sample_gwas_df, sample_eqtl_df):
        """Bokeh backend should handle eQTL panel with effect sizes without error."""
        plotter = LocusZoomPlotter(species=None, backend="bokeh", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE1",
        )

        assert fig is not None

    def test_bokeh_eqtl_without_effects(self, sample_gwas_df, sample_eqtl_df_no_effect):
        """Bokeh backend should handle eQTL panel without effect sizes."""
        plotter = LocusZoomPlotter(species=None, backend="bokeh", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df_no_effect,
            eqtl_gene="GENE1",
        )

        assert fig is not None

    def test_bokeh_finemapping(self, sample_gwas_df, sample_finemapping_df):
        """Bokeh backend should handle fine-mapping panel without error."""
        plotter = LocusZoomPlotter(species=None, backend="bokeh", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
        )

        assert fig is not None

    def test_plotly_combined_eqtl_finemapping(
        self, sample_gwas_df, sample_eqtl_df, sample_finemapping_df
    ):
        """Plotly backend should handle both eQTL and fine-mapping panels together."""
        plotter = LocusZoomPlotter(species=None, backend="plotly", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE1",
            finemapping_df=sample_finemapping_df,
        )

        assert fig is not None

    def test_bokeh_combined_eqtl_finemapping(
        self, sample_gwas_df, sample_eqtl_df, sample_finemapping_df
    ):
        """Bokeh backend should handle both eQTL and fine-mapping panels together."""
        plotter = LocusZoomPlotter(species=None, backend="bokeh", log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE1",
            finemapping_df=sample_finemapping_df,
        )

        assert fig is not None

    def test_eqtl_chr_filtering(self, sample_gwas_df):
        """Test that eQTL panel filters by chromosome, not just position."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        # Create eQTL data with chr column, some on wrong chromosome
        eqtl_df = pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],  # All in region 1e6-2e6
                "p_value": [1e-6, 1e-4, 0.01],
                "gene": ["GENE1", "GENE1", "GENE1"],
                "effect_size": [0.5, -0.3, 0.8],
                "chr": ["1", "2", "1"],  # Middle one is on chr2
            }
        )

        # Plot for chr 1 - should only include 2 eQTLs (not the chr2 one)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=eqtl_df,
            eqtl_gene="GENE1",
        )

        assert fig is not None
        plt.close(fig)

    def test_eqtl_gene_without_gene_column_no_gene_in_label(self, sample_gwas_df):
        """Test that eqtl_gene without 'gene' column doesn't label as gene-filtered.

        Bug fix: Previously, label would show "eQTL (GENE1)" even when filtering
        didn't occur because the DataFrame lacked a 'gene' column.

        Warning is logged to stderr (see "Captured stderr call" in test output).
        loguru doesn't integrate with pytest's caplog/capsys fixtures directly.
        """
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        # Create eQTL data WITHOUT gene column
        eqtl_df_no_gene_col = pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],
                "p_value": [1e-6, 1e-4, 0.01],
                # No "gene" column - so filtering can't occur
            }
        )

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
            eqtl_df=eqtl_df_no_gene_col,
            eqtl_gene="GENE1",  # Specified but can't filter
        )

        # Verify the eQTL panel axes don't have "(GENE1)" in any labels
        # The panel is axes[1] (after GWAS panel at axes[0])
        axes = fig.get_axes()
        eqtl_ax = axes[1]  # eQTL panel

        # Check that no legend entry contains "(GENE1)"
        legend = eqtl_ax.get_legend()
        if legend:
            for text in legend.get_texts():
                assert "(GENE1)" not in text.get_text(), (
                    f"Label incorrectly shows gene filter: {text.get_text()}"
                )

        assert fig is not None
        plt.close(fig)


class TestPheWASPlot:
    """Tests for plot_phewas method."""

    def test_plot_phewas_basic(self):
        """Test basic PheWAS plot generation."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        phewas_df = pd.DataFrame(
            {
                "phenotype": ["Height", "BMI", "T2D", "CAD", "HDL"],
                "p_value": [1e-15, 0.05, 1e-8, 1e-3, 1e-10],
                "category": [
                    "Anthropometric",
                    "Anthropometric",
                    "Metabolic",
                    "Cardiovascular",
                    "Metabolic",
                ],
            }
        )

        fig = plotter.plot_phewas(phewas_df, variant_id="rs12345")

        assert fig is not None
        plt.close(fig)

    def test_plot_phewas_with_effect(self):
        """Test PheWAS plot with effect sizes."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        phewas_df = pd.DataFrame(
            {
                "phenotype": ["Height", "BMI", "T2D"],
                "p_value": [1e-15, 0.05, 1e-8],
                "category": ["Anthropometric", "Anthropometric", "Metabolic"],
                "effect_size": [0.5, -0.1, 0.3],
            }
        )

        fig = plotter.plot_phewas(
            phewas_df, variant_id="rs12345", effect_col="effect_size"
        )

        assert fig is not None
        plt.close(fig)

    def test_plot_phewas_no_category(self):
        """Test PheWAS plot without category column."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        phewas_df = pd.DataFrame(
            {
                "phenotype": ["Height", "BMI", "T2D"],
                "p_value": [1e-15, 0.05, 1e-8],
            }
        )

        fig = plotter.plot_phewas(phewas_df, variant_id="rs12345")

        assert fig is not None
        plt.close(fig)


class TestForestPlot:
    """Tests for plot_forest method."""

    def test_plot_forest_basic(self):
        """Test basic forest plot generation."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        forest_df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Study C", "Meta-analysis"],
                "effect": [0.5, 0.3, 0.6, 0.45],
                "ci_lower": [0.2, 0.0, 0.3, 0.35],
                "ci_upper": [0.8, 0.6, 0.9, 0.55],
            }
        )

        fig = plotter.plot_forest(forest_df, variant_id="rs12345")

        assert fig is not None
        plt.close(fig)

    def test_plot_forest_with_weights(self):
        """Test forest plot with study weights."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        forest_df = pd.DataFrame(
            {
                "study": ["Study A", "Study B", "Meta"],
                "effect": [0.5, 0.3, 0.4],
                "ci_lower": [0.2, 0.0, 0.3],
                "ci_upper": [0.8, 0.6, 0.5],
                "weight": [30, 70, 100],
            }
        )

        fig = plotter.plot_forest(forest_df, variant_id="rs12345", weight_col="weight")

        assert fig is not None
        plt.close(fig)

    def test_plot_forest_custom_null_value(self):
        """Test forest plot with custom null value (e.g., OR=1)."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        forest_df = pd.DataFrame(
            {
                "study": ["Study A", "Study B"],
                "effect": [1.5, 0.9],
                "ci_lower": [1.1, 0.6],
                "ci_upper": [1.9, 1.3],
            }
        )

        fig = plotter.plot_forest(
            forest_df,
            variant_id="rs12345",
            null_value=1.0,
            effect_label="Odds Ratio",
        )

        assert fig is not None
        plt.close(fig)


class TestRecombinationDownloadErrors:
    """Tests for recombination map error handling.

    These tests verify that when recombination maps are unavailable,
    the plotter gracefully handles None return values and allows
    plotting to continue without recombination overlay.

    Note: Detailed error handling (network, I/O, OS errors) is tested
    in test_recombination.py at the ensure_recomb_maps level.
    """

    @pytest.fixture
    def plotter(self):
        """Create a plotter instance for testing download errors."""
        return LocusZoomPlotter(species="canine", log_level="DEBUG")

    def test_ensure_recomb_maps_returns_none_propagates(self, plotter):
        """When ensure_recomb_maps returns None, plotter._ensure_recomb_maps returns None."""
        with patch(
            "pylocuszoom.plotter.ensure_recomb_maps", return_value=None
        ) as mock_ensure:
            result = plotter._ensure_recomb_maps()
            assert result is None
            mock_ensure.assert_called_once_with(species="canine", data_dir=None)

    def test_plotting_continues_without_recomb_maps(self, plotter):
        """Plotting should succeed even when recombination maps are unavailable."""
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom.plotter.ensure_recomb_maps", return_value=None):
            # Should not raise, just skip recombination overlay
            fig = plotter.plot(
                gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=True,  # Requested but unavailable
            )
            assert fig is not None
            plt.close(fig)


class TestPvalueTransformation:
    """Tests for p-value transformation helper."""

    def test_transform_pvalues_adds_neglog10p_column(self):
        """Helper creates neglog10p column from p-values."""
        df = pd.DataFrame({"pval": [0.01, 0.001, 1e-8]})
        plotter = LocusZoomPlotter()

        result = plotter._transform_pvalues(df.copy(), "pval")

        assert "neglog10p" in result.columns
        assert result["neglog10p"].iloc[0] == pytest.approx(2.0)  # -log10(0.01)
        assert result["neglog10p"].iloc[1] == pytest.approx(3.0)  # -log10(0.001)
        assert result["neglog10p"].iloc[2] == pytest.approx(8.0)  # -log10(1e-8)

    def test_transform_pvalues_clips_extreme_values(self):
        """Extremely small p-values are clipped to avoid -inf."""
        df = pd.DataFrame({"pval": [1e-350, 0.0]})  # Would be -inf without clipping
        plotter = LocusZoomPlotter()

        result = plotter._transform_pvalues(df.copy(), "pval")

        # Should be clipped to 1e-300, giving ~300
        assert result["neglog10p"].iloc[0] == pytest.approx(300.0)
        assert result["neglog10p"].iloc[1] == pytest.approx(300.0)
        assert not np.isinf(result["neglog10p"]).any()


class TestPlotterDelegation:
    """Tests for plotter delegation to specialized classes."""

    @pytest.fixture
    def plotter(self):
        """Create a LocusZoomPlotter instance."""
        return LocusZoomPlotter(species="canine")

    def test_manhattan_delegation_preserves_species(self, plotter):
        """Test that species is passed to ManhattanPlotter."""
        assert plotter._manhattan_plotter.species == "canine"

    def test_stats_delegation_preserves_threshold(self, plotter):
        """Test that threshold is passed to StatsPlotter."""
        assert (
            plotter._stats_plotter.genomewide_threshold == plotter.genomewide_threshold
        )

    def test_plotter_delegates_to_plot_finemapping(self):
        """Test that plot_stacked delegates finemapping to module function."""
        with patch("pylocuszoom.plotter.plot_finemapping") as mock_plot_fm:
            plotter = LocusZoomPlotter(
                species="canine", backend="matplotlib", log_level=None
            )

            gwas_df = pd.DataFrame(
                {
                    "ps": [1000, 2000],
                    "p_wald": [0.01, 0.001],
                }
            )
            fm_df = pd.DataFrame(
                {
                    "pos": [1000, 2000],
                    "pip": [0.5, 0.3],
                    "cs": [1, 1],
                }
            )

            with patch.object(plotter, "_get_recomb_for_region", return_value=None):
                try:
                    plotter.plot_stacked(
                        [gwas_df],
                        chrom=1,
                        start=0,
                        end=3000,
                        finemapping_df=fm_df,
                    )
                except Exception:
                    pass  # We just want to verify the call was made

            mock_plot_fm.assert_called()


def test_plotter_uses_ensure_recomb_maps():
    """Test that LocusZoomPlotter._ensure_recomb_maps delegates to module function."""
    with patch(
        "pylocuszoom.plotter.ensure_recomb_maps",
        return_value=Path("/mock/recomb"),
    ) as mock_ensure:
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        # Call the method that should delegate to ensure_recomb_maps
        result = plotter._ensure_recomb_maps()

        mock_ensure.assert_called_once_with(species="canine", data_dir=None)
        assert result == Path("/mock/recomb")


class TestDeprecationWarnings:
    """Tests for deprecation warnings on wrapper methods."""

    def test_plot_manhattan_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_manhattan emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame(
            {
                "chrom": [1, 1],
                "pos": [1000, 2000],
                "p": [0.01, 0.001],
            }
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_manhattan(df)
            except Exception:
                pass  # May fail on rendering, we just want the warning

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "ManhattanPlotter" in str(deprecation_warnings[0].message)

    def test_plot_qq_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_qq emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame({"p": [0.01, 0.001, 0.5]})

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_qq(df)
            except Exception:
                pass

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "ManhattanPlotter" in str(deprecation_warnings[0].message)

    def test_plot_manhattan_stacked_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_manhattan_stacked emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame(
            {
                "chrom": [1, 1],
                "pos": [1000, 2000],
                "p": [0.01, 0.001],
            }
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_manhattan_stacked([df])
            except Exception:
                pass

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "ManhattanPlotter" in str(deprecation_warnings[0].message)

    def test_plot_manhattan_qq_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_manhattan_qq emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame(
            {
                "chrom": [1, 1],
                "pos": [1000, 2000],
                "p": [0.01, 0.001],
            }
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_manhattan_qq(df)
            except Exception:
                pass

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "ManhattanPlotter" in str(deprecation_warnings[0].message)

    def test_plot_manhattan_qq_stacked_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_manhattan_qq_stacked emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame(
            {
                "chrom": [1, 1],
                "pos": [1000, 2000],
                "p": [0.01, 0.001],
            }
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_manhattan_qq_stacked([df])
            except Exception:
                pass

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "ManhattanPlotter" in str(deprecation_warnings[0].message)

    def test_plot_phewas_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_phewas emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame(
            {
                "phenotype": ["A", "B"],
                "p_value": [0.01, 0.001],
                "category": ["cat1", "cat1"],
            }
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_phewas(df, variant_id="rs123")
            except Exception:
                pass

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "StatsPlotter" in str(deprecation_warnings[0].message)

    def test_plot_forest_warns_deprecated(self):
        """Test that LocusZoomPlotter.plot_forest emits deprecation warning."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        df = pd.DataFrame(
            {
                "study": ["A", "B"],
                "effect": [0.5, 0.3],
                "ci_lower": [0.2, 0.1],
                "ci_upper": [0.8, 0.5],
            }
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                plotter.plot_forest(df, variant_id="rs123")
            except Exception:
                pass

            deprecation_warnings = [
                x for x in w if issubclass(x.category, DeprecationWarning)
            ]
            assert len(deprecation_warnings) >= 1
            assert "StatsPlotter" in str(deprecation_warnings[0].message)

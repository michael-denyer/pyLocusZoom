"""Tests for LocusZoomPlotter class."""

from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import settings as hyp_settings

from pylocuszoom._plotter_utils import transform_pvalues
from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend
from pylocuszoom.backends.plotly_backend import PlotlyBackend
from pylocuszoom.plotter import LocusZoomPlotter
from tests.strategies import gwas_dataframes


class TestBackendIntegration:
    """Tests for backend protocol integration."""

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
        """plot() with backend='plotly' produces a plotly Figure.

        This also implicitly confirms plot() routes through the backend
        protocol: if plot() bypassed the backend and called matplotlib
        directly, the returned object would be a matplotlib Figure and
        this isinstance check would fail.
        """
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

    def test_matplotlib_plot_renders_expected_artists(self, sample_gwas_df):
        """plot() renders scatter points, a significance line, and axis labels.

        Assertions query the rendered axes directly (scatter collections,
        line objects, y-label text) rather than mocking backend methods
        and counting calls — per CLAUDE.md's observable-outputs rule.
        """
        plotter = LocusZoomPlotter(species="canine")

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        try:
            ax = fig.axes[0]

            # At least one scatter collection (association points)
            assert len(ax.collections) >= 1, (
                "expected at least one scatter collection on association axes"
            )

            # At least one line (the significance threshold)
            lines = ax.get_lines()
            assert len(lines) >= 1, "expected a significance threshold line"

            # Y-axis labelled with -log10(p)
            ylabel = ax.get_ylabel()
            assert "log" in ylabel.lower() or "-log" in ylabel.lower(), (
                f"expected -log10(p) y-label, got {ylabel!r}"
            )

            # X-limits cover the requested region
            xlim = ax.get_xlim()
            assert xlim[0] <= 1_000_000 and xlim[1] >= 2_000_000
        finally:
            plt.close(fig)

    def test_plot_stacked_renders_two_panels(self, sample_gwas_df):
        """plot_stacked() with two GWAS inputs produces two association panels."""
        plotter = LocusZoomPlotter(species="canine")

        fig = plotter.plot_stacked(
            [sample_gwas_df, sample_gwas_df.copy()],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        try:
            # Each GWAS gets its own association axes; gene track axes
            # may also be present but at minimum we need two scatter axes.
            scatter_axes = [ax for ax in fig.axes if ax.collections]
            assert len(scatter_axes) >= 2, (
                f"expected >=2 axes with scatter data, got {len(scatter_axes)}"
            )
        finally:
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
        rng = np.random.default_rng(42)
        n_snps = 50
        positions = np.sort(rng.integers(1000000, 2000000, n_snps))
        return pd.DataFrame(
            {
                "rs": [f"rs{i}" for i in range(n_snps)],
                "chr": [1] * n_snps,
                "ps": positions,
                "p_wald": rng.uniform(1e-10, 1, n_snps),
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
        """Empty GWAS DataFrame should raise ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        empty_df = pd.DataFrame(columns=["rs", "chr", "ps", "p_wald"])
        with pytest.raises(ValidationError, match="empty"):
            plotter.plot(
                empty_df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

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
        df["R2"] = np.random.default_rng(0).uniform(0, 1, len(df))

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

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
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
                show_recombination=False,
            )

            mock_ld.assert_called_once()
            plt.close(fig)

    def test_empty_ld_output_is_downgraded_to_warning(self, plotter, caplog):
        """An empty-output PlinkError (singleton lead SNP) should not abort.

        Singleton lead SNPs with no LD neighbours in the window are a real
        scenario; plotter.plot() catches only this specific PlinkError and
        continues without LD colouring, leaving a warning in the log.
        """
        from pylocuszoom.exceptions import EmptyLDOutputError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.side_effect = EmptyLDOutputError(
                "PLINK produced an empty LD output for lead SNP 'rs1'."
            )
            fig = plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_pos=1100000,
                ld_reference_file="/path/to/genotypes",
                show_recombination=False,
            )
        assert fig is not None
        plt.close(fig)

    def test_stacked_plot_downgrades_empty_ld_output(self, plotter):
        """Single and stacked plots share the recoverable LD policy."""
        from pylocuszoom.exceptions import EmptyLDOutputError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.side_effect = EmptyLDOutputError("no LD neighbours")
            fig = plotter.plot_stacked(
                [df],
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=[1100000],
                ld_reference_files=["/path/to/genotypes"],
                show_recombination=False,
            )

        assert fig is not None
        plt.close(fig)

    def test_plink_misconfiguration_propagates_through_plot(self, plotter):
        """A non-empty-output PlinkError must surface — it means PLINK is broken.

        Regression boundary: the catch in plotter.plot() is narrow on purpose.
        Timeout, non-zero exit, and "output file missing after success" all
        indicate real misconfiguration and should reach the caller.
        """
        from pylocuszoom.exceptions import PlinkError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.side_effect = PlinkError(
                "PLINK LD calculation failed (exit code 2): bad bfile"
            )
            with pytest.raises(PlinkError, match="exit code"):
                plotter.plot(
                    df,
                    chrom=1,
                    start=1000000,
                    end=2000000,
                    lead_pos=1100000,
                    ld_reference_file="/path/to/genotypes",
                )


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

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
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

    def test_transform_pvalues_filters_invalid_range(self, plotter):
        """_transform_pvalues filters out-of-range p-values (< 0 or > 1)."""
        import io

        from loguru import logger as loguru_logger

        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1700000, 1900000],
                "p_wald": [0.001, 0.5, -0.1, 1.5],
            }
        )

        # Capture log output
        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="WARNING",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        try:
            result = transform_pvalues(df.copy(), "p_wald")
        finally:
            loguru_logger.remove(handler_id)

        # Only valid rows (0.001, 0.5) should remain
        assert len(result) == 2
        assert "neglog10p" in result.columns
        assert np.isfinite(result["neglog10p"]).all()
        assert (result["neglog10p"] > 0).all()

        # Warning logged about 2 invalid values
        log_output = log_capture.getvalue()
        assert "2 p-values outside [0, 1]" in log_output

    def test_transform_pvalues_filters_nan(self, plotter):
        """_transform_pvalues filters NaN p-values."""
        import io

        from loguru import logger as loguru_logger

        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [0.001, np.nan, 0.5],
            }
        )

        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="WARNING",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        try:
            result = transform_pvalues(df.copy(), "p_wald")
        finally:
            loguru_logger.remove(handler_id)

        # NaN row filtered out
        assert len(result) == 2
        assert not result["p_wald"].isna().any()

        # Warning logged
        log_output = log_capture.getvalue()
        assert "1 NaN p-values" in log_output

    def test_transform_pvalues_clips_very_small(self, plotter):
        """_transform_pvalues clips very small p-values to 1e-300."""
        import io

        from loguru import logger as loguru_logger

        df = pd.DataFrame(
            {
                "ps": [1100000],
                "p_wald": [1e-310],  # Smaller than 1e-300
            }
        )

        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="DEBUG",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        try:
            result = transform_pvalues(df.copy(), "p_wald")
        finally:
            loguru_logger.remove(handler_id)

        # Should be clipped to 1e-300, giving -log10(1e-300) = 300
        assert len(result) == 1
        assert result["neglog10p"].iloc[0] == pytest.approx(300.0)
        assert not np.isinf(result["neglog10p"]).any()

        # Debug log about clipping
        log_output = log_capture.getvalue()
        assert "Clipping" in log_output

    def test_transform_pvalues_preserves_valid_data(self, plotter):
        """_transform_pvalues passes through all-valid data unchanged."""
        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [0.001, 0.5, 1e-8],
            }
        )

        result = transform_pvalues(df.copy(), "p_wald")

        assert len(result) == 3
        assert result["neglog10p"].iloc[0] == pytest.approx(3.0)
        assert result["neglog10p"].iloc[1] == pytest.approx(-np.log10(0.5))
        assert result["neglog10p"].iloc[2] == pytest.approx(8.0)

    def test_plot_stacked_lead_detection_excludes_out_of_range(self, plotter):
        """Lead SNP auto-detection should exclude out-of-range p-values.

        Regression test: lead detection only filtered NaN but not out-of-range
        p-values, so a lead SNP with p=-0.1 could be selected and then removed
        by _transform_pvalues, causing missing lead highlighting.
        """
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                # rs1 has smallest absolute value but is invalid (negative)
                # rs3 should be selected as lead (smallest valid p-value)
                "p_wald": [-0.1, 0.5, 0.001],
            }
        )

        # Should not raise — lead detection should skip the invalid p-value
        fig = plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        plt.close(fig)

    def test_plot_stacked_all_invalid_pvalues(self, plotter):
        """plot_stacked should handle region with all out-of-range p-values."""
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "chr": [1, 1],
                "ps": [1100000, 1500000],
                "p_wald": [-0.1, 1.5],  # All invalid
            }
        )

        fig = plotter.plot_stacked(
            [gwas_df],
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

        result = transform_pvalues(df.copy(), "pval")

        assert "neglog10p" in result.columns
        assert result["neglog10p"].iloc[0] == pytest.approx(2.0)  # -log10(0.01)
        assert result["neglog10p"].iloc[1] == pytest.approx(3.0)  # -log10(0.001)
        assert result["neglog10p"].iloc[2] == pytest.approx(8.0)  # -log10(1e-8)

    def test_transform_pvalues_clips_extreme_values(self):
        """Extremely small p-values are clipped to avoid -inf."""
        df = pd.DataFrame({"pval": [1e-350, 0.0]})  # Would be -inf without clipping

        result = transform_pvalues(df.copy(), "pval")

        # Should be clipped to 1e-300, giving ~300
        assert result["neglog10p"].iloc[0] == pytest.approx(300.0)
        assert result["neglog10p"].iloc[1] == pytest.approx(300.0)
        assert not np.isinf(result["neglog10p"]).any()


class TestPlotterDelegation:
    """Tests for plotter delegation to specialized classes."""

    def test_composer_delegates_to_plot_finemapping(self):
        """The composer forwards prepared fine-mapping data to its renderer.

        Pins the dispatch contract at the layer that now owns panel rendering.
        Asserts on the DataFrame handed to plot_finemapping, not just that it
        was called.
        """
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

        with (
            patch.object(plotter, "_get_recomb_for_region", return_value=None),
            patch("pylocuszoom._regional.plot_finemapping") as mock_plot_fm,
        ):
            fig = plotter.plot_stacked(
                [gwas_df],
                chrom=1,
                start=0,
                end=3000,
                finemapping_df=fm_df,
            )
            try:
                mock_plot_fm.assert_called_once()
                # Signature: plot_finemapping(backend, ax, df, ...). The
                # DataFrame is positional index 2 (or kwarg 'df').
                args, kwargs = mock_plot_fm.call_args
                forwarded = kwargs.get("df")
                if forwarded is None and len(args) >= 3:
                    forwarded = args[2]
                assert forwarded is not None, (
                    "plot_finemapping was called without the finemapping DataFrame"
                )
                assert list(forwarded["pos"]) == [1000, 2000]
                assert list(forwarded["pip"]) == [0.5, 0.3]
            finally:
                plt.close(fig)


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


# =============================================================================
# Property-Based Tests (Hypothesis)
# =============================================================================


class TestPlotterProperties:
    """Property-based tests for plotter crash-resistance."""

    @hyp_settings(max_examples=15, deadline=None)
    @given(gwas_dataframes(min_snps=3, max_snps=50))
    def test_plot_renders_valid_data_matplotlib(self, df):
        """plot() should render any valid GWAS data without crashing (matplotlib)."""
        plotter = LocusZoomPlotter(species="canine")
        chrom = df["chr"].iloc[0]
        start = int(df["ps"].min())
        end = int(df["ps"].max())

        fig = plotter.plot(
            df,
            chrom=chrom,
            start=start,
            end=end,
            show_recombination=False,
        )

        assert fig is not None
        plt.close(fig)

    @hyp_settings(max_examples=10, deadline=None)
    @given(gwas_dataframes(min_snps=3, max_snps=30))
    def test_plot_renders_valid_data_plotly(self, df):
        """plot() should render any valid GWAS data without crashing (plotly)."""
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(species="canine", backend="plotly")
        chrom = df["chr"].iloc[0]
        start = int(df["ps"].min())
        end = int(df["ps"].max())

        fig = plotter.plot(
            df,
            chrom=chrom,
            start=start,
            end=end,
            show_recombination=False,
        )

        assert isinstance(fig, go.Figure)

    @hyp_settings(max_examples=10, deadline=None)
    @given(gwas_dataframes(min_snps=3, max_snps=30))
    def test_plot_renders_valid_data_bokeh(self, df):
        """plot() should render any valid GWAS data without crashing (bokeh)."""
        plotter = LocusZoomPlotter(species="canine", backend="bokeh")
        chrom = df["chr"].iloc[0]
        start = int(df["ps"].min())
        end = int(df["ps"].max())

        fig = plotter.plot(
            df,
            chrom=chrom,
            start=start,
            end=end,
            show_recombination=False,
        )

        assert fig is not None


class TestPlotStackedProperties:
    """Property-based tests for stacked plots."""

    @hyp_settings(max_examples=10, deadline=None)
    @given(gwas_dataframes(min_snps=5, max_snps=30))
    def test_plot_stacked_with_duplicated_data(self, df):
        """plot_stacked() should handle multiple identical DataFrames."""
        plotter = LocusZoomPlotter(species="canine")
        chrom = df["chr"].iloc[0]
        start = int(df["ps"].min())
        end = int(df["ps"].max())

        fig = plotter.plot_stacked(
            [df, df.copy()],
            chrom=chrom,
            start=start,
            end=end,
            show_recombination=False,
        )

        assert fig is not None
        plt.close(fig)


# =============================================================================
# LD Heatmap Integration Tests
# =============================================================================


class TestLDHeatmapIntegration:
    """Tests for LD heatmap integration in LocusZoomPlotter."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS DataFrame with positions matching LD heatmap SNPs."""
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "chr": [1, 1, 1, 1, 1],
                "ps": [1000000, 1000500, 1001000, 1001500, 1002000],
                "p_wald": [1e-8, 1e-6, 1e-4, 1e-3, 0.05],
            }
        )

    @pytest.fixture
    def sample_ld_heatmap_data(self):
        """Sample LD heatmap data matching calculate_pairwise_ld return format.

        Returns:
            Tuple of (DataFrame, list[str]): LD matrix and SNP IDs.
        """
        ld_matrix = pd.DataFrame(
            np.array(
                [
                    [1.0, 0.9, 0.7, 0.4, 0.2],
                    [0.9, 1.0, 0.8, 0.5, 0.3],
                    [0.7, 0.8, 1.0, 0.6, 0.4],
                    [0.4, 0.5, 0.6, 1.0, 0.7],
                    [0.2, 0.3, 0.4, 0.7, 1.0],
                ]
            ),
            index=["rs1", "rs2", "rs3", "rs4", "rs5"],
            columns=["rs1", "rs2", "rs3", "rs4", "rs5"],
        )
        snp_ids = ["rs1", "rs2", "rs3", "rs4", "rs5"]
        return ld_matrix, snp_ids

    @pytest.fixture
    def sample_genes_df(self):
        """Sample genes for testing stacked plots with gene track."""
        return pd.DataFrame(
            {
                "chr": ["1", "1"],
                "start": [1000200, 1001200],
                "end": [1000800, 1001800],
                "gene_name": ["GENE1", "GENE2"],
                "strand": ["+", "-"],
            }
        )

    def test_plot_with_ld_heatmap_renders_two_panels(
        self, sample_gwas_df, sample_ld_heatmap_data
    ):
        """When ld_heatmap_df and ld_heatmap_snp_ids provided, figure has heatmap panel."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should have at least 2 axes (association + heatmap)
        assert len(fig.axes) >= 2
        plt.close(fig)

    def test_plot_with_ld_heatmap_aligns_x_coordinates(
        self, sample_gwas_df, sample_ld_heatmap_data
    ):
        """Heatmap SNPs render at their genomic positions from GWAS data."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Verify the heatmap panel exists and has correct x-axis range
        # The heatmap x-axis should span the genomic positions of SNPs
        axes = fig.axes
        assert len(axes) >= 2
        # Heatmap is axes[1] (after association at axes[0], colorbar may be axes[2])
        heatmap_ax = axes[1]
        xlim = heatmap_ax.get_xlim()
        # X-axis should be in genomic coordinate range
        assert xlim[0] < 1003000, f"Heatmap xlim[0]={xlim[0]} should be < 1003000"
        assert xlim[1] > 999000, f"Heatmap xlim[1]={xlim[1]} should be > 999000"
        plt.close(fig)

    def test_plot_stacked_with_ld_heatmap_at_bottom(
        self, sample_gwas_df, sample_ld_heatmap_data, sample_genes_df
    ):
        """In stacked plots, heatmap appears below gene track (at very bottom)."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            genes_df=sample_genes_df,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should have 3 panels: GWAS, gene track, heatmap
        axes = fig.axes
        assert len(axes) >= 3
        plt.close(fig)

    def test_ld_heatmap_filters_to_region(self, sample_gwas_df, sample_ld_heatmap_data):
        """SNPs outside [start, end] are filtered from heatmap."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Narrow region that only includes some SNPs
        # rs1 at 1000000, rs2 at 1000500 - both in [1000000, 1001000)
        # rs3 at 1001000 - at boundary
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1000000,
            end=1001000,  # Only includes rs1, rs2, rs3
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error even with partial overlap
        assert fig is not None
        plt.close(fig)

    def test_ld_heatmap_empty_overlap_logs_warning(
        self, sample_gwas_df, sample_ld_heatmap_data, capfd
    ):
        """When no SNPs in heatmap overlap with region, warning logged and no heatmap panel added."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level="WARNING")

        # Region that doesn't overlap with any GWAS SNP positions
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=5000000,  # Far outside GWAS positions
            end=6000000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error - heatmap simply not added
        assert fig is not None
        plt.close(fig)

    def test_ld_heatmap_height_parameter(self, sample_gwas_df, sample_ld_heatmap_data):
        """ld_heatmap_height parameter controls panel height ratio."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Test with different height values
        fig1 = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
            ld_heatmap_height=0.1,  # Smaller
        )

        fig2 = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
            ld_heatmap_height=0.5,  # Larger
        )

        # Both should render successfully
        assert fig1 is not None
        assert fig2 is not None
        plt.close(fig1)
        plt.close(fig2)

    def test_ld_heatmap_lead_snp_highlight(
        self, sample_gwas_df, sample_ld_heatmap_data
    ):
        """Lead SNP is highlighted in heatmap when lead_pos specified and SNP found."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # rs1 is at position 1000000
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            lead_pos=1000000,  # rs1
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should render with lead SNP highlight (visual check)
        assert fig is not None
        plt.close(fig)

    # Backend-specific tests

    def test_ld_heatmap_matplotlib_backend(
        self, sample_gwas_df, sample_ld_heatmap_data
    ):
        """Verify matplotlib figure has correct panel count and axes."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should have at least 2 main axes (association + heatmap)
        # Plus possible colorbar axis
        assert len(fig.axes) >= 2
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_ld_heatmap_plotly_backend(self, sample_gwas_df, sample_ld_heatmap_data):
        """Verify plotly figure has heatmap trace at correct row."""
        import plotly.graph_objects as go

        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend="plotly", log_level=None)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        assert isinstance(fig, go.Figure)
        # Check that figure has data traces
        assert len(fig.data) > 0

    def test_ld_heatmap_bokeh_backend(self, sample_gwas_df, sample_ld_heatmap_data):
        """Verify bokeh layout contains heatmap."""
        from bokeh.models.layouts import Column

        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend="bokeh", log_level=None)

        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Bokeh returns Column layout
        assert isinstance(fig, Column)
        # Should have multiple children (panels)
        assert len(fig.children) >= 2

    # Edge case tests

    def test_ld_heatmap_single_snp_in_region(self, sample_gwas_df):
        """Test with single SNP in filtered region (graceful handling)."""
        # Create a matrix with 5 SNPs but only 1 will be in region
        ld_matrix = pd.DataFrame(
            np.array(
                [
                    [1.0, 0.9, 0.7, 0.4, 0.2],
                    [0.9, 1.0, 0.8, 0.5, 0.3],
                    [0.7, 0.8, 1.0, 0.6, 0.4],
                    [0.4, 0.5, 0.6, 1.0, 0.7],
                    [0.2, 0.3, 0.4, 0.7, 1.0],
                ]
            ),
            index=["rs1", "rs2", "rs3", "rs4", "rs5"],
            columns=["rs1", "rs2", "rs3", "rs4", "rs5"],
        )
        snp_ids = ["rs1", "rs2", "rs3", "rs4", "rs5"]

        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Very narrow region that only includes rs1 (at position 1000000)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=999999,
            end=1000001,  # Only rs1 at 1000000
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error
        assert fig is not None
        plt.close(fig)

    def test_ld_heatmap_lead_snp_not_in_heatmap(
        self, sample_gwas_df, sample_ld_heatmap_data
    ):
        """Test lead SNP not in heatmap SNPs (no highlighting, no error)."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Create GWAS with a lead SNP not in the heatmap
        gwas_with_extra = sample_gwas_df.copy()
        gwas_with_extra = pd.concat(
            [
                gwas_with_extra,
                pd.DataFrame(
                    {
                        "rs": ["rs_extra"],
                        "chr": [1],
                        "ps": [1000100],  # Different position
                        "p_wald": [1e-10],  # Most significant
                    }
                ),
            ],
            ignore_index=True,
        )

        fig = plotter.plot(
            gwas_with_extra,
            chrom=1,
            start=999000,
            end=1003000,
            lead_pos=1000100,  # rs_extra - not in heatmap
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error
        assert fig is not None
        plt.close(fig)

    def test_ld_heatmap_missing_snp_ids_raises_error(
        self, sample_gwas_df, sample_ld_heatmap_data
    ):
        """Test that providing ld_heatmap_df without ld_heatmap_snp_ids raises error."""
        ld_matrix, _ = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with pytest.raises(ValueError, match="ld_heatmap_snp_ids is required"):
            plotter.plot(
                sample_gwas_df,
                chrom=1,
                start=999000,
                end=1003000,
                show_recombination=False,
                ld_heatmap_df=ld_matrix,
                # ld_heatmap_snp_ids not provided
            )


class TestHighlightHeatmapSnpBackendProtocol:
    """Tests that highlight_heatmap_snp works on all backends."""

    def test_matplotlib_highlight_heatmap_snp(self):
        """Matplotlib backend highlight_heatmap_snp creates rectangle patches."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(6, 4)
        )
        ax = axes[0]

        # Should not raise
        backend.highlight_heatmap_snp(ax, fig, snp_idx=2, n_snps=5)

        # Verify patches were added (3 for row + 2 for column = 5 total)
        rect_patches = [
            p for p in ax.patches if hasattr(p, "get_edgecolor") and not p.get_fill()
        ]
        assert len(rect_patches) == 5
        plt.close(fig)

    def test_plotly_highlight_heatmap_snp(self):
        """Plotly backend highlight_heatmap_snp adds shapes to figure."""
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, panels = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(6, 4)
        )
        ax = panels[0]

        # Should not raise
        backend.highlight_heatmap_snp(ax, fig, snp_idx=2, n_snps=5)

        # Verify shapes were added
        shapes = fig.layout.shapes
        assert len(shapes) == 5  # 3 for row + 2 for column

    def test_bokeh_highlight_heatmap_snp(self):
        """Bokeh backend highlight_heatmap_snp adds rect glyphs to figure."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        layout, figures = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(6, 4)
        )
        ax = figures[0]
        renderers_before = len(ax.renderers)

        # Should not raise
        backend.highlight_heatmap_snp(ax, layout, snp_idx=2, n_snps=5)

        # Verify renderer was added (batched into single rect call)
        assert len(ax.renderers) == renderers_before + 1

        # Verify the batched renderer has correct number of cells
        # snp_idx=2, n_snps=5: row cells (0,2),(1,2),(2,2) + col cells (2,3),(2,4) = 5
        renderer = ax.renderers[-1]
        assert len(renderer.data_source.data["x"]) == 5


class TestStackedPlotLeadDetectionCrossChrom:
    """Regression: plot_stacked() lead auto-detection must filter by chromosome.

    Pre-fix bug: when lead_positions is None, the auto-detect loop did
    `df[(df[pos_col] >= start) & (df[pos_col] <= end)]` with no chrom
    filter. On a multi-chromosome GWAS DataFrame, the strongest p-value
    in [start, end] could come from a different chromosome than the
    plotted region, anchoring the diamond marker to the wrong locus.
    """

    @pytest.fixture
    def plotter(self):
        return LocusZoomPlotter(species="canine", log_level=None)

    def test_lead_autodetect_filters_by_chrom(self, plotter):
        """Lead position must come from the requested chromosome."""
        # Two chromosomes share a position range. The strongest p-value
        # is on chr2 at position 1_500_000, but we are plotting chr1.
        # Pre-fix: lead = chr2's strongest hit. Post-fix: lead = chr1's.
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4"],
                "chr": [1, 1, 2, 2],
                "ps": [1_200_000, 1_800_000, 1_500_000, 1_900_000],
                "p_wald": [1e-5, 1e-3, 1e-12, 1e-10],
            }
        )

        # Patch the lead-aware code to capture what got computed.
        captured = {}
        composer = plotter._regional_composer
        original = composer.render_association_scatter

        def spy(ax, df, pos_col, ld_col, lead_pos, *args, **kwargs):
            captured.setdefault("lead_positions", []).append(lead_pos)
            return original(ax, df, pos_col, ld_col, lead_pos, *args, **kwargs)

        composer.render_association_scatter = spy
        try:
            plotter.plot_stacked(
                [gwas_df],
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                pos_col="ps",
                p_col="p_wald",
                show_recombination=False,
            )
        finally:
            composer.render_association_scatter = original

        assert captured["lead_positions"] == [1_200_000], (
            "Lead must be chr1's strongest hit (1_200_000), not chr2's (1_500_000)"
        )


class TestLeadPosBoundary:
    """Pins the lead_pos=1 boundary.

    Smallest valid position (``ge=1`` in config) reaches
    ``render_association_scatter`` intact, and the public API rejects ``0``
    (1-based genomic coords).
    """

    def test_lead_pos_one_reaches_plot_association(self):
        """lead_pos=1 (smallest valid position) reaches the association scatter."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs_lead", "rs2", "rs3"],
                "ps": [1, 100_000, 200_000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        captured = {}
        composer = plotter._regional_composer
        original = composer.render_association_scatter

        def spy(ax, df, pos_col, ld_col, lead_pos, *args, **kwargs):
            captured["lead_pos"] = lead_pos
            return original(ax, df, pos_col, ld_col, lead_pos, *args, **kwargs)

        composer.render_association_scatter = spy
        try:
            plotter.plot(
                gwas_df,
                chrom=1,
                start=1,
                end=300_000,
                pos_col="ps",
                p_col="p_wald",
                lead_pos=1,
                show_recombination=False,
            )
        finally:
            composer.render_association_scatter = original

        assert captured["lead_pos"] == 1, (
            "lead_pos=1 must pass through; falsy-check regression would drop it to None"
        )

    def test_lead_pos_zero_rejected_at_api(self):
        """Public API enforces genomic coords are 1-based; lead_pos=0 rejected."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "ps": [100_000, 200_000],
                "p_wald": [1e-8, 1e-5],
            }
        )

        from pydantic import ValidationError as PydanticValidationError

        with pytest.raises(PydanticValidationError, match="greater than or equal to 1"):
            plotter.plot(
                gwas_df,
                chrom=1,
                start=1,
                end=300_000,
                pos_col="ps",
                p_col="p_wald",
                lead_pos=0,
                show_recombination=False,
            )

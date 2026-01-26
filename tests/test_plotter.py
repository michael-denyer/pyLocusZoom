"""Tests for LocusZoomPlotter class."""

from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.plotter import LocusZoomPlotter


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
                eqtl_df=bad_eqtl_df,
                show_recombination=False,
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

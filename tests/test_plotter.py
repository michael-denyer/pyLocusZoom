"""Tests for LocusZoomPlotter class."""

from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.plotter import LocusZoomPlotter


class TestLocusZoomPlotterInit:
    """Tests for LocusZoomPlotter initialization."""

    def test_default_species_is_dog(self):
        """Default species should be dog."""
        plotter = LocusZoomPlotter()
        assert plotter.species == "dog"

    def test_custom_species(self):
        """Should accept custom species."""
        plotter = LocusZoomPlotter(species="cat")
        assert plotter.species == "cat"

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
        return LocusZoomPlotter(species="dog")

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
        return LocusZoomPlotter(species="dog", plink_path="/mock/plink")

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

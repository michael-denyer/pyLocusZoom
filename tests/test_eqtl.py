"""Tests for eQTL module filtering and colocalization functions."""

import pandas as pd
import pytest

from pylocuszoom.eqtl import (
    EQTLValidationError,
    calculate_colocalization_overlap,
    filter_eqtl_by_gene,
    filter_eqtl_by_region,
    get_eqtl_genes,
    prepare_eqtl_for_plotting,
    validate_eqtl_df,
)


class TestValidateEqtlDf:
    """Tests for validate_eqtl_df function."""

    def test_valid_eqtl_passes(self):
        """Valid eQTL DataFrame passes validation."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
            }
        )
        # Should not raise
        validate_eqtl_df(df)

    def test_missing_position_column_fails(self):
        """Missing position column raises EQTLValidationError."""
        df = pd.DataFrame(
            {
                "p_value": [1e-6, 0.01],
            }
        )
        with pytest.raises(EQTLValidationError):
            validate_eqtl_df(df)

    def test_missing_pvalue_column_fails(self):
        """Missing p-value column raises EQTLValidationError."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
            }
        )
        with pytest.raises(EQTLValidationError):
            validate_eqtl_df(df)

    def test_custom_column_names(self):
        """Custom column names work correctly."""
        df = pd.DataFrame(
            {
                "position": [1000000, 1001000],
                "pval": [1e-6, 0.01],
            }
        )
        validate_eqtl_df(df, pos_col="position", p_col="pval")


class TestFilterEqtlByGene:
    """Tests for filter_eqtl_by_gene function."""

    @pytest.fixture
    def eqtl_df(self):
        """Sample eQTL data with multiple genes."""
        return pd.DataFrame(
            {
                "pos": [1000000, 1001000, 1002000, 1003000],
                "p_value": [1e-6, 0.01, 1e-8, 0.05],
                "gene": ["BRCA1", "BRCA1", "TP53", "TP53"],
            }
        )

    def test_filters_to_single_gene(self, eqtl_df):
        """Filters to single gene correctly."""
        result = filter_eqtl_by_gene(eqtl_df, gene="BRCA1")
        assert len(result) == 2
        assert all(result["gene"] == "BRCA1")

    def test_returns_copy(self, eqtl_df):
        """Returns a copy, not a view."""
        result = filter_eqtl_by_gene(eqtl_df, gene="BRCA1")
        result.loc[result.index[0], "p_value"] = 0.99
        assert eqtl_df["p_value"].iloc[0] == 1e-6

    def test_missing_gene_column_raises_error(self):
        """Missing gene column raises EQTLValidationError."""
        df = pd.DataFrame(
            {
                "pos": [1000000],
                "p_value": [1e-6],
            }
        )
        with pytest.raises(EQTLValidationError, match="gene"):
            filter_eqtl_by_gene(df, gene="BRCA1")

    def test_custom_gene_column(self):
        """Custom gene column name works."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
                "gene_name": ["BRCA1", "TP53"],
            }
        )
        result = filter_eqtl_by_gene(df, gene="BRCA1", gene_col="gene_name")
        assert len(result) == 1

    def test_no_matching_gene_returns_empty(self, eqtl_df):
        """No matching gene returns empty DataFrame."""
        result = filter_eqtl_by_gene(eqtl_df, gene="NONEXISTENT")
        assert len(result) == 0
        assert isinstance(result, pd.DataFrame)


class TestFilterEqtlByRegion:
    """Tests for filter_eqtl_by_region function."""

    @pytest.fixture
    def eqtl_df(self):
        """Sample eQTL data spanning a region."""
        return pd.DataFrame(
            {
                "chr": [1, 1, 1, 2],
                "pos": [1000000, 1500000, 2000000, 1500000],
                "p_value": [1e-6, 0.01, 1e-8, 0.05],
            }
        )

    def test_filters_to_region(self, eqtl_df):
        """Filters to genomic region correctly."""
        result = filter_eqtl_by_region(eqtl_df, chrom=1, start=1200000, end=1800000)
        assert len(result) == 1
        assert result["pos"].iloc[0] == 1500000

    def test_filters_by_chromosome(self, eqtl_df):
        """Filters by chromosome when chr column present."""
        result = filter_eqtl_by_region(eqtl_df, chrom=1, start=0, end=3000000)
        assert len(result) == 3
        assert all(result["chr"] == 1)

    def test_no_chr_column_filters_position_only(self):
        """Works without chromosome column when chrom_col is empty string."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "p_value": [1e-6, 0.01, 1e-8],
            }
        )
        result = filter_eqtl_by_region(
            df, chrom=1, start=1200000, end=1800000, chrom_col=""
        )
        assert len(result) == 1


class TestGetEqtlGenes:
    """Tests for get_eqtl_genes function."""

    def test_returns_unique_sorted_genes(self):
        """Returns unique genes sorted alphabetically."""
        df = pd.DataFrame(
            {
                "gene": ["TP53", "BRCA1", "TP53", "EGFR"],
            }
        )
        result = get_eqtl_genes(df)
        assert result == ["BRCA1", "EGFR", "TP53"]

    def test_missing_gene_column_returns_empty_list(self):
        """Missing gene column returns empty list."""
        df = pd.DataFrame(
            {
                "pos": [1000000],
            }
        )
        result = get_eqtl_genes(df)
        assert result == []

    def test_handles_nan_values(self):
        """NaN values are excluded."""
        df = pd.DataFrame(
            {
                "gene": ["BRCA1", None, "TP53"],
            }
        )
        result = get_eqtl_genes(df)
        assert result == ["BRCA1", "TP53"]


class TestCalculateColocalizationOverlap:
    """Tests for calculate_colocalization_overlap function."""

    @pytest.fixture
    def gwas_df(self):
        """Sample GWAS data."""
        return pd.DataFrame(
            {
                "ps": [1000000, 1001000, 1002000, 1003000],
                "p_wald": [1e-6, 0.01, 1e-8, 0.05],
            }
        )

    @pytest.fixture
    def eqtl_df(self):
        """Sample eQTL data with some overlap."""
        return pd.DataFrame(
            {
                "pos": [1000000, 1002000, 1004000],  # 1000000 and 1002000 overlap
                "p_value": [1e-7, 1e-6, 1e-8],
            }
        )

    def test_finds_overlapping_significant_snps(self, gwas_df, eqtl_df):
        """Finds SNPs significant in both datasets."""
        result = calculate_colocalization_overlap(gwas_df, eqtl_df, p_threshold=1e-5)
        # Only 1000000 and 1002000 are significant in both
        assert len(result) == 2
        assert 1000000 in result["ps"].values
        assert 1002000 in result["ps"].values

    def test_no_overlap_returns_empty(self, gwas_df):
        """No overlapping positions returns empty DataFrame."""
        eqtl_no_overlap = pd.DataFrame(
            {
                "pos": [9000000, 9001000],
                "p_value": [1e-10, 1e-10],
            }
        )
        result = calculate_colocalization_overlap(
            gwas_df, eqtl_no_overlap, p_threshold=1e-5
        )
        assert len(result) == 0

    def test_custom_column_names(self):
        """Custom column names work correctly."""
        gwas = pd.DataFrame(
            {
                "position": [1000000],
                "pvalue": [1e-10],
            }
        )
        eqtl = pd.DataFrame(
            {
                "bp": [1000000],
                "p": [1e-10],
            }
        )
        result = calculate_colocalization_overlap(
            gwas,
            eqtl,
            gwas_pos_col="position",
            eqtl_pos_col="bp",
            gwas_p_col="pvalue",
            eqtl_p_col="p",
        )
        assert len(result) == 1


class TestPrepareEqtlForPlotting:
    """Tests for prepare_eqtl_for_plotting function."""

    def test_adds_neglog10p_column(self):
        """Adds neglog10p column to DataFrame."""
        df = pd.DataFrame(
            {
                "pos": [1000000],
                "p_value": [1e-6],
            }
        )
        result = prepare_eqtl_for_plotting(df)
        assert "neglog10p" in result.columns
        assert result["neglog10p"].iloc[0] == pytest.approx(6.0)

    def test_filters_by_gene_when_specified(self):
        """Filters by gene when gene parameter provided."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
                "gene": ["BRCA1", "TP53"],
            }
        )
        result = prepare_eqtl_for_plotting(df, gene="BRCA1")
        assert len(result) == 1

    def test_filters_by_region_when_specified(self):
        """Filters by region when all region params provided."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "p_value": [1e-6, 0.01, 1e-8],
            }
        )
        result = prepare_eqtl_for_plotting(df, chrom=1, start=1200000, end=1800000)
        assert len(result) == 1

    def test_returns_copy_not_modifying_original(self):
        """Returns a copy without modifying original DataFrame."""
        df = pd.DataFrame(
            {
                "pos": [1000000],
                "p_value": [1e-6],
            }
        )
        result = prepare_eqtl_for_plotting(df)
        assert "neglog10p" not in df.columns
        assert "neglog10p" in result.columns

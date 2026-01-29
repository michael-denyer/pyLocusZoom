"""Tests for SNP label placement module."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.labels import add_snp_labels


class TestAddSnpLabels:
    """Tests for add_snp_labels function."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results with neglog10p calculated."""
        np.random.seed(42)
        return pd.DataFrame(
            {
                "rs": [
                    "rs1",
                    "rs2",
                    "rs3",
                    "rs4",
                    "rs5",
                    "rs6",
                    "rs7",
                    "rs8",
                    "rs9",
                    "rs10",
                ],
                "ps": [
                    1100000,
                    1200000,
                    1300000,
                    1400000,
                    1500000,
                    1600000,
                    1700000,
                    1800000,
                    1900000,
                    2000000,
                ],
                "p_wald": [1e-8, 1e-5, 1e-3, 1e-6, 1e-9, 1e-4, 1e-7, 1e-2, 1e-10, 1e-1],
                "neglog10p": [8, 5, 3, 6, 9, 4, 7, 2, 10, 1],
            }
        )

    @pytest.fixture
    def sample_genes_df(self):
        """Sample gene annotations."""
        return pd.DataFrame(
            {
                "chr": [1, 1, 1],
                "start": [1050000, 1450000, 1850000],
                "end": [1150000, 1550000, 1950000],
                "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
            }
        )

    def test_adds_labels_to_plot(self, sample_gwas_df):
        """Should add text labels to the plot."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=3)

        assert len(texts) == 3
        plt.close(fig)

    def test_labels_top_n_snps(self, sample_gwas_df):
        """Should label only the top N most significant SNPs."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=5)

        assert len(texts) == 5
        # Top 5 by neglog10p should be: rs9 (10), rs5 (9), rs1 (8), rs7 (7), rs4 (6)
        plt.close(fig)

    def test_uses_snp_id_by_default(self, sample_gwas_df):
        """Should use SNP ID (rs number) as default label."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=1)

        # Top SNP is rs9 with neglog10p=10
        assert "rs9" in texts[0].get_text()
        plt.close(fig)

    def test_ignores_genes_df_uses_snp_id(self, sample_gwas_df, sample_genes_df):
        """Should always use SNP ID, even when genes_df provided."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        # rs9 is at position 1900000 - should still show rs9, not gene name
        texts = add_snp_labels(
            ax, sample_gwas_df, label_top_n=1, genes_df=sample_genes_df, chrom=1
        )

        assert "rs9" in texts[0].get_text()
        plt.close(fig)

    def test_truncates_long_labels(self, sample_gwas_df):
        """Should truncate labels longer than max_label_length."""
        # Create SNP with very long name
        df = sample_gwas_df.copy()
        df.loc[df["neglog10p"] == 10, "rs"] = "rs_very_long_identifier_name"

        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(ax, df, label_top_n=1, max_label_length=15)

        label_text = texts[0].get_text()
        assert len(label_text) <= 15
        assert label_text.endswith("...")
        plt.close(fig)

    def test_raises_without_neglog10p_column(self, sample_gwas_df):
        """Should raise error if neglog10p column is missing."""
        df = sample_gwas_df.drop(columns=["neglog10p"])

        fig, ax = plt.subplots()

        with pytest.raises(ValueError, match="neglog10p"):
            add_snp_labels(ax, df, label_top_n=1)

        plt.close(fig)

    def test_handles_custom_column_names(self):
        """Should work with custom column names."""
        df = pd.DataFrame(
            {
                "snp_id": ["var1", "var2", "var3"],
                "position": [1000, 2000, 3000],
                "log_pval": [5, 10, 3],
            }
        )

        fig, ax = plt.subplots()
        ax.scatter(df["position"], df["log_pval"])

        texts = add_snp_labels(
            ax,
            df,
            pos_col="position",
            neglog10p_col="log_pval",
            rs_col="snp_id",
            label_top_n=2,
        )

        assert len(texts) == 2
        # Top should be var2 with log_pval=10
        assert "var2" in texts[0].get_text()
        plt.close(fig)

    def test_returns_empty_list_for_zero_labels(self, sample_gwas_df):
        """Should return empty list when label_top_n=0."""
        fig, ax = plt.subplots()

        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=0)

        assert len(texts) == 0
        plt.close(fig)


class TestAdjustTextWarning:
    """Test warning is logged when adjustText is unavailable."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results with neglog10p calculated."""
        np.random.seed(42)
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "ps": [1100000, 1200000, 1300000, 1400000, 1500000],
                "neglog10p": [8, 5, 3, 6, 9],
            }
        )

    def test_warning_logged_when_adjusttext_unavailable(self, sample_gwas_df):
        """Verify warning is logged when adjustText import fails."""
        import builtins
        import io
        from unittest.mock import patch

        from loguru import logger as loguru_logger

        from pylocuszoom.logging import logger

        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        # Capture log output with a custom sink
        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="WARNING",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        # Store the original __import__
        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "adjustText":
                raise ImportError("No module named 'adjustText'")
            return original_import(name, *args, **kwargs)

        # Ensure logging is enabled and patch __import__ to fail for adjustText
        logger.enable("WARNING")
        try:
            with patch.object(builtins, "__import__", side_effect=mock_import):
                # Need multiple labels to trigger the adjustText code path
                texts = add_snp_labels(
                    ax=ax,
                    df=sample_gwas_df,
                    label_top_n=3,
                )
        finally:
            loguru_logger.remove(handler_id)
            plt.close(fig)

        # Should still return labels even without adjustText
        assert len(texts) == 3

        # Check that warning was logged about adjustText
        log_output = log_capture.getvalue()
        assert "adjustText" in log_output, (
            f"Expected warning about adjustText, got: {log_output}"
        )

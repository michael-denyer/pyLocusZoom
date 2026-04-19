"""Tests for SNP label placement module."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.labels import add_snp_labels, adjust_snp_labels


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

    def test_filters_labels_near_lead_snp(self):
        """Non-lead SNPs within min_label_distance of lead are excluded."""
        # Region is 1Mb wide, 5% threshold = 50kb
        # Lead at 1500000, nearby SNP at 1510000 (10kb away = within 5%)
        # Distant SNP at 1800000 (300kb away = outside 5%)
        df = pd.DataFrame(
            {
                "rs": ["lead", "nearby", "distant"],
                "ps": [1500000, 1510000, 1800000],
                "neglog10p": [20, 15, 10],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(
            ax,
            df,
            label_top_n=3,
            lead_pos=1500000,
            region_span=1000000,
        )

        labels = [t.get_text() for t in texts]
        assert "lead" in labels
        assert "distant" in labels
        assert "nearby" not in labels
        plt.close(fig)

    def test_no_filtering_without_lead_pos(self):
        """All top N labels shown when lead_pos is not provided."""
        df = pd.DataFrame(
            {
                "rs": ["a", "b", "c"],
                "ps": [1500000, 1510000, 1800000],
                "neglog10p": [20, 15, 10],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(ax, df, label_top_n=3)

        assert len(texts) == 3
        plt.close(fig)

    def test_lead_snp_always_labeled(self):
        """Lead SNP within top-N is always labeled even when neighbors filtered."""
        df = pd.DataFrame(
            {
                "rs": ["lead", "neighbor1", "neighbor2"],
                "ps": [1500000, 1500100, 1499900],
                "neglog10p": [20, 18, 16],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(
            ax,
            df,
            label_top_n=3,
            lead_pos=1500000,
            region_span=1000000,
        )

        labels = [t.get_text() for t in texts]
        assert "lead" in labels
        assert len(labels) == 1  # Only lead survives
        plt.close(fig)

    def test_snp_exactly_at_threshold_survives(self):
        """SNP exactly at min_label_distance threshold is kept (>= boundary)."""
        # 5% of 1Mb = 50kb; SNP at exactly 50kb from lead should survive
        df = pd.DataFrame(
            {
                "rs": ["lead", "boundary"],
                "ps": [1500000, 1550000],
                "neglog10p": [20, 15],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(
            ax, df, label_top_n=2, lead_pos=1500000, region_span=1000000
        )

        labels = [t.get_text() for t in texts]
        assert "boundary" in labels
        plt.close(fig)

    def test_no_filtering_when_region_span_zero(self):
        """region_span=0 disables filtering and logs a warning."""
        df = pd.DataFrame(
            {
                "rs": ["a", "b", "c"],
                "ps": [1500000, 1510000, 1800000],
                "neglog10p": [20, 15, 10],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(ax, df, label_top_n=3, lead_pos=1500000, region_span=0)

        assert len(texts) == 3
        plt.close(fig)

    def test_warns_when_lead_pos_without_region_span(self, caplog):
        """Warning logged when lead_pos set but region_span missing."""
        df = pd.DataFrame(
            {
                "rs": ["a", "b"],
                "ps": [1500000, 1800000],
                "neglog10p": [20, 10],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        with caplog.at_level("WARNING"):
            texts = add_snp_labels(
                ax, df, label_top_n=2, lead_pos=1500000, region_span=None
            )

        # All labels kept (filtering disabled)
        assert len(texts) == 2
        plt.close(fig)

    def test_custom_min_label_distance(self):
        """Custom min_label_distance widens the exclusion zone."""
        # 10% of 1Mb = 100kb; SNP 60kb away should be filtered at 10%
        # but would survive at the default 5% (50kb threshold)
        df = pd.DataFrame(
            {
                "rs": ["lead", "mid_range", "far"],
                "ps": [1500000, 1560000, 1800000],
                "neglog10p": [20, 15, 10],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        texts = add_snp_labels(
            ax,
            df,
            label_top_n=3,
            lead_pos=1500000,
            region_span=1000000,
            min_label_distance=0.10,
        )

        labels = [t.get_text() for t in texts]
        assert "lead" in labels
        assert "far" in labels
        assert "mid_range" not in labels
        plt.close(fig)

    def test_invalid_min_label_distance_raises(self):
        """min_label_distance outside [0, 1] raises ValueError."""
        df = pd.DataFrame(
            {
                "rs": ["a", "b"],
                "ps": [1500000, 1800000],
                "neglog10p": [20, 10],
            }
        )
        fig, ax = plt.subplots()
        ax.scatter(df["ps"], df["neglog10p"])

        with pytest.raises(ValueError, match="min_label_distance"):
            add_snp_labels(
                ax,
                df,
                label_top_n=2,
                lead_pos=1500000,
                region_span=1000000,
                min_label_distance=-0.5,
            )
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


class TestDeferredAdjustment:
    """Test deferred label adjustment for proper axis limits handling."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Sample GWAS results with neglog10p calculated."""
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "ps": [1100000, 1200000, 1300000, 1400000, 1500000],
                "neglog10p": [8, 5, 3, 6, 9],
            }
        )

    def test_adjust_false_skips_adjustment(self, sample_gwas_df):
        """When adjust=False, labels are added but not adjusted."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=3, adjust=False)

        assert len(texts) == 3
        plt.close(fig)

    def test_adjust_snp_labels_can_be_called_separately(self, sample_gwas_df):
        """adjust_snp_labels can be called after setting axis limits."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        # Add labels without adjustment
        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=3, adjust=False)

        # Set axis limits
        ax.set_xlim(1000000, 1600000)
        ax.set_ylim(0, 12)

        # Now adjust labels
        adjust_snp_labels(ax, texts)

        assert len(texts) == 3
        plt.close(fig)

    def test_adjust_snp_labels_handles_empty_list(self):
        """adjust_snp_labels handles empty text list gracefully."""
        fig, ax = plt.subplots()

        # Should not raise
        adjust_snp_labels(ax, [])

        plt.close(fig)

    def test_adjust_snp_labels_handles_single_label(self, sample_gwas_df):
        """adjust_snp_labels handles single label (skips adjustText)."""
        fig, ax = plt.subplots()
        ax.scatter(sample_gwas_df["ps"], sample_gwas_df["neglog10p"])

        texts = add_snp_labels(ax, sample_gwas_df, label_top_n=1, adjust=False)

        # Should not raise - adjustText is skipped for single label
        adjust_snp_labels(ax, texts)

        assert len(texts) == 1
        plt.close(fig)


class TestNearLeadFilterBackfill:
    """Regression: filter near-lead SNPs BEFORE nlargest, not after.

    Pre-fix bug: nlargest(label_top_n) was taken first, then the
    proximity mask dropped near-lead neighbors. On a strong peak the
    top N are clustered around the lead, so the user got 1 label
    instead of label_top_n.
    """

    def test_near_lead_filter_backfills_to_label_top_n(self):
        """All top-N labels should fit in budget after near-lead filter."""
        # Lead at 1500000; many strong-but-near neighbors that should be
        # excluded; weaker far-away SNPs that should backfill the budget.
        df = pd.DataFrame(
            {
                "rs": [f"rs{i}" for i in range(10)],
                "ps": [
                    1500000,  # lead
                    1500100,  # near (excluded)
                    1500200,  # near (excluded)
                    1500300,  # near (excluded)
                    1500400,  # near (excluded)
                    1100000,  # far - should be picked
                    1200000,  # far - should be picked
                    1300000,  # far - should be picked
                    1700000,  # far - should be picked
                    1900000,  # far - should be picked (5th non-lead)
                ],
                "neglog10p": [
                    20.0,  # lead - strongest
                    19.0,
                    18.0,
                    17.0,
                    16.0,
                    5.0,
                    4.0,
                    3.0,
                    2.5,
                    2.0,
                ],
            }
        )

        fig, ax = plt.subplots()
        try:
            texts = add_snp_labels(
                ax,
                df,
                label_top_n=5,
                lead_pos=1500000,
                region_span=1_000_000,
                min_label_distance=0.05,
                adjust=False,
            )
            # Pre-fix: only the lead survives the post-filter mask -> 1 label.
            # Post-fix: lead + 4 farthest backfill to 5.
            assert len(texts) == 5
            label_strs = [t.get_text() for t in texts]
            assert "rs0" in label_strs  # lead always present
        finally:
            plt.close(fig)

    def test_near_lead_filter_preserves_lead_when_only_neighbors_compete(self):
        """If only near-lead points exist, lead is still labeled."""
        df = pd.DataFrame(
            {
                "rs": ["lead", "near1", "near2"],
                "ps": [1500000, 1500100, 1500200],
                "neglog10p": [20.0, 19.0, 18.0],
            }
        )
        fig, ax = plt.subplots()
        try:
            texts = add_snp_labels(
                ax,
                df,
                label_top_n=3,
                lead_pos=1500000,
                region_span=1_000_000,
                adjust=False,
            )
            label_strs = [t.get_text() for t in texts]
            assert "lead" in label_strs
            # Near neighbors are excluded; only lead is eligible.
            assert len(texts) == 1
        finally:
            plt.close(fig)

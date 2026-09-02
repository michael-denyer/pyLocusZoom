"""Tests for the LD enrichment policy in _ld_plotting.py.

Every branch here decides whether to reach for PLINK at all. The branches that
decline return the caller's frame untouched, so none of these tests needs a
PLINK binary or a subprocess.
"""

import pandas as pd

from pylocuszoom._ld_plotting import enrich_with_ld

ARGS = {
    "pos_col": "pos",
    "rs_col": "rs",
    "start": 1_000_000,
    "end": 2_000_000,
    "plink_path": None,
    "species": "canine",
}


class TestEnrichWithLDDeclines:
    """The cases that return the frame unchanged, without calling PLINK."""

    def test_no_reference_file_returns_the_frame_unchanged(self, tiny_regional_gwas_df):
        """Without a reference panel there is nothing to compute LD against."""
        result, ld_col = enrich_with_ld(
            tiny_regional_gwas_df,
            reference_file=None,
            lead_pos=1_500_000,
            ld_col=None,
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col is None

    def test_no_lead_position_returns_the_frame_unchanged(self, tiny_regional_gwas_df):
        """LD is measured against a lead SNP, so no lead means no LD."""
        result, ld_col = enrich_with_ld(
            tiny_regional_gwas_df,
            reference_file="/nonexistent/panel",
            lead_pos=None,
            ld_col=None,
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col is None

    def test_existing_ld_column_wins_over_the_reference_panel(
        self, tiny_regional_gwas_df
    ):
        """A caller-supplied LD column is kept and PLINK is not consulted."""
        result, ld_col = enrich_with_ld(
            tiny_regional_gwas_df,
            reference_file="/nonexistent/panel",
            lead_pos=1_500_000,
            ld_col="R2",
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col == "R2"

    def test_missing_rs_column_warns_and_returns_the_frame(
        self, tiny_regional_gwas_df, warning_records
    ):
        """PLINK addresses variants by ID, so no ID column means no LD."""
        without_ids = tiny_regional_gwas_df.drop(columns=["rs"])

        result, ld_col = enrich_with_ld(
            without_ids,
            reference_file="/nonexistent/panel",
            lead_pos=1_500_000,
            ld_col=None,
            **ARGS,
        )

        assert result is without_ids
        assert ld_col is None
        assert any("'rs' not found" in record for record in warning_records)

    def test_lead_position_absent_from_the_frame_warns(
        self, tiny_regional_gwas_df, warning_records
    ):
        """A lead position with no matching row cannot name a lead SNP."""
        result, ld_col = enrich_with_ld(
            tiny_regional_gwas_df,
            reference_file="/nonexistent/panel",
            lead_pos=1_234_567,
            ld_col=None,
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col is None
        assert any("1234567 not found" in record for record in warning_records)


class TestEnrichWithLDMerges:
    """The success path merges PLINK's R2 column onto the caller's frame."""

    def test_r2_is_merged_onto_matching_variant_ids(
        self, monkeypatch, tiny_regional_gwas_df
    ):
        """Each variant gains the R2 that PLINK reported for its ID."""
        import pylocuszoom._ld_plotting as module

        monkeypatch.setattr(
            module,
            "calculate_ld",
            lambda **kwargs: pd.DataFrame(
                {"SNP": ["rs1", "rs2", "rs3"], "R2": [0.9, 1.0, 0.2]}
            ),
        )

        result, ld_col = enrich_with_ld(
            tiny_regional_gwas_df,
            reference_file="/panel",
            lead_pos=1_500_000,
            ld_col=None,
            **ARGS,
        )

        assert ld_col == "R2"
        assert list(result["R2"]) == [0.9, 1.0, 0.2]

    def test_empty_ld_output_warns_and_keeps_the_frame(
        self, monkeypatch, tiny_regional_gwas_df, warning_records
    ):
        """A region PLINK finds no pairs in degrades to an uncoloured plot."""
        import pylocuszoom._ld_plotting as module
        from pylocuszoom.exceptions import EmptyLDOutputError

        def raise_empty(**kwargs):
            raise EmptyLDOutputError("no pairs in window")

        monkeypatch.setattr(module, "calculate_ld", raise_empty)

        result, ld_col = enrich_with_ld(
            tiny_regional_gwas_df,
            reference_file="/panel",
            lead_pos=1_500_000,
            ld_col=None,
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col is None
        assert any("no pairs in window" in record for record in warning_records)

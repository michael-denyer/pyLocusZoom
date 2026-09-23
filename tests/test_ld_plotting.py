"""Tests for the LD enrichment policy in _ld_plotting.py.

Every branch here decides whether to reach for PLINK at all. The branches that
decline return the caller's frame untouched, so none of these tests needs a
PLINK binary or a subprocess.
"""

import pandas as pd
import pytest

from pylocuszoom._ld_plotting import enrich_with_ld

ARGS = {
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
            lead_index=1,
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
            lead_index=None,
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
            lead_index=1,
            ld_col="R2",
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col == "R2"


class TestEnrichWithLDLookup:
    """The success path assigns PLINK's R2 values to the selected rows."""

    def test_r2_is_assigned_to_matching_variant_ids(
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
            lead_index=1,
            ld_col=None,
            **ARGS,
        )

        assert ld_col == "R2"
        assert list(result["R2"]) == [0.9, 1.0, 0.2]

    def test_lookup_preserves_rows_and_replaces_old_r2(self, monkeypatch):
        frame = pd.DataFrame(
            {"rs": ["b", "a", "missing", "a"], "R2": [-1.0] * 4},
            index=[8, 3, 9, 5],
        )
        original = frame.copy(deep=True)
        monkeypatch.setattr(
            "pylocuszoom._ld_plotting.calculate_ld",
            lambda **kwargs: pd.DataFrame({"SNP": ["a", "b"], "R2": [0.5, 1.0]}),
        )
        result, ld_col = enrich_with_ld(
            frame, reference_file="/panel", lead_index=8, ld_col=None, **ARGS
        )
        expected = original.assign(R2=[1.0, 0.5, float("nan"), 0.5])
        assert ld_col == "R2"
        pd.testing.assert_frame_equal(result, expected)
        pd.testing.assert_frame_equal(frame, original)

    def test_duplicate_reference_ids_are_rejected(
        self, monkeypatch, tiny_regional_gwas_df
    ):
        monkeypatch.setattr(
            "pylocuszoom._ld_plotting.calculate_ld",
            lambda **kwargs: pd.DataFrame({"SNP": ["rs1", "rs1"], "R2": [0.5, 1.0]}),
        )
        with pytest.raises(ValueError):
            enrich_with_ld(
                tiny_regional_gwas_df,
                reference_file="/panel",
                lead_index=0,
                ld_col=None,
                **ARGS,
            )

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
            lead_index=1,
            ld_col=None,
            **ARGS,
        )

        assert result is tiny_regional_gwas_df
        assert ld_col is None
        assert any("no pairs in window" in record for record in warning_records)

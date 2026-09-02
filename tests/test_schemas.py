"""Tests for the load-time contracts in schemas.py.

Tests edge cases and error paths for GWAS, eQTL, fine-mapping, and gene
annotation DataFrame validation, and file path validation.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import LoaderValidationError
from pylocuszoom.schemas import Family, Tier, spec
from pylocuszoom.validation import check

# =============================================================================
# GWAS validation edge cases
# =============================================================================


class TestGWASValidationEdgeCases:
    """Edge cases for the GWAS load-time contract beyond happy/error paths."""

    def test_valid_gwas_passes(self):
        """A minimal GWAS frame with an rs column satisfies the contract."""
        df = pd.DataFrame({"ps": [100], "p_wald": [0.5], "rs": ["rs1"]})
        check(df, spec(Family.GWAS, Tier.LOAD))

    def test_custom_column_names_validated(self):
        """Custom pos_col and p_col names are respected during validation."""
        df = pd.DataFrame({"position": [100], "pvalue": [0.5]})
        check(df, spec(Family.GWAS, Tier.LOAD, pos_col="position", p_col="pvalue"))

    def test_custom_column_names_missing_raises(self):
        """Missing custom column names raise LoaderValidationError."""
        df = pd.DataFrame({"ps": [100], "p_wald": [0.5]})
        with pytest.raises(LoaderValidationError, match=r"Missing columns.*'pos'"):
            check(df, spec(Family.GWAS, Tier.LOAD, pos_col="pos", p_col="p_wald"))

    def test_both_columns_missing_reports_both(self):
        """Error message mentions both missing columns when both are absent."""
        df = pd.DataFrame({"chr": [1], "rs": ["rs1"]})
        with pytest.raises(LoaderValidationError, match="ps") as exc_info:
            check(df, spec(Family.GWAS, Tier.LOAD))
        assert "p_wald" in str(exc_info.value)

    def test_zero_position_fails(self):
        """Position of exactly zero is rejected (non-positive)."""
        df = pd.DataFrame({"ps": [0], "p_wald": [0.5]})
        with pytest.raises(LoaderValidationError, match=r"'ps': 1 values <= 0"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_pvalue_exactly_one_passes(self):
        """P-value of exactly 1.0 is valid (range is (0, 1])."""
        df = pd.DataFrame({"ps": [100], "p_wald": [1.0]})
        check(df, spec(Family.GWAS, Tier.LOAD))

    def test_nan_pvalue_fails(self):
        """NaN p-values are caught."""
        df = pd.DataFrame({"ps": [100, 200], "p_wald": [0.05, float("nan")]})
        with pytest.raises(LoaderValidationError, match=r"'p_wald' has 1 null values"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_multiple_numeric_errors_reported(self):
        """Multiple validation errors accumulated and reported at once."""
        df = pd.DataFrame(
            {
                "ps": ["chr1:100"],
                "p_wald": ["not_a_number"],
            }  # non-numeric  # non-numeric
        )
        with pytest.raises(LoaderValidationError, match="must be numeric") as exc_info:
            check(df, spec(Family.GWAS, Tier.LOAD))
        msg = str(exc_info.value)
        assert "ps" in msg
        assert "p_wald" in msg

    def test_single_row_valid(self):
        """Single-row DataFrame passes validation."""
        df = pd.DataFrame({"ps": [42], "p_wald": [1e-10]})
        check(df, spec(Family.GWAS, Tier.LOAD))

    def test_extra_columns_ignored(self):
        """Extra columns beyond required ones don't cause failures."""
        df = pd.DataFrame(
            {"ps": [100], "p_wald": [0.5], "rs": ["rs1"], "beta": [0.3], "se": [0.1]}
        )
        check(df, spec(Family.GWAS, Tier.LOAD))


# =============================================================================
# eQTL validation edge cases
# =============================================================================


class TestEQTLValidationEdgeCases:
    """Edge cases for the eQTL load-time contract."""

    def test_valid_eqtl_passes(self):
        """A minimal eQTL frame satisfies the contract."""
        df = pd.DataFrame({"pos": [100], "p_value": [0.05], "gene": ["BRCA1"]})
        check(df, spec(Family.EQTL, Tier.LOAD))

    def test_non_numeric_pos_fails(self):
        """Non-numeric position column is rejected."""
        df = pd.DataFrame({"pos": ["chr1:100"], "p_value": [0.05], "gene": ["BRCA1"]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.EQTL, Tier.LOAD))

    def test_non_numeric_pvalue_fails(self):
        """Non-numeric p_value column is rejected."""
        df = pd.DataFrame({"pos": [100], "p_value": ["significant"], "gene": ["BRCA1"]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.EQTL, Tier.LOAD))

    def test_zero_pvalue_fails(self):
        """P-value of 0 is rejected (outside (0, 1])."""
        df = pd.DataFrame({"pos": [100], "p_value": [0.0], "gene": ["BRCA1"]})
        with pytest.raises(LoaderValidationError, match=r"'p_value': 1 values <= 0"):
            check(df, spec(Family.EQTL, Tier.LOAD))

    def test_missing_multiple_columns(self):
        """All missing required columns listed in error."""
        df = pd.DataFrame({"extra": [1]})
        with pytest.raises(LoaderValidationError) as exc_info:
            check(df, spec(Family.EQTL, Tier.LOAD))
        msg = str(exc_info.value)
        assert "pos" in msg
        assert "p_value" in msg
        assert "gene" in msg

    def test_non_positive_pos_fails(self):
        """Non-positive positions are rejected."""
        df = pd.DataFrame(
            {"pos": [-5, 100], "p_value": [0.05, 0.1], "gene": ["A", "B"]}
        )
        with pytest.raises(LoaderValidationError, match=r"'pos': 1 values <= 0"):
            check(df, spec(Family.EQTL, Tier.LOAD))


# =============================================================================
# Fine-mapping validation edge cases
# =============================================================================


class TestFinemappingValidationEdgeCases:
    """Edge cases for the fine-mapping load-time contract."""

    def test_pip_zero_and_one_valid(self):
        """PIP of 0.0 and 1.0 are both valid (range [0, 1])."""
        df = pd.DataFrame({"pos": [100, 200], "pip": [0.0, 1.0]})
        check(df, spec(Family.FINEMAPPING, Tier.LOAD))

    def test_non_numeric_pip_fails(self):
        """Non-numeric pip column is rejected."""
        df = pd.DataFrame({"pos": [100], "pip": ["high"]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.FINEMAPPING, Tier.LOAD))

    def test_non_numeric_pos_fails(self):
        """Non-numeric pos column is rejected."""
        df = pd.DataFrame({"pos": ["abc"], "pip": [0.5]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.FINEMAPPING, Tier.LOAD))

    def test_missing_both_columns(self):
        """Both pos and pip missing reported."""
        df = pd.DataFrame({"other": [1]})
        with pytest.raises(LoaderValidationError) as exc_info:
            check(df, spec(Family.FINEMAPPING, Tier.LOAD))
        msg = str(exc_info.value)
        assert "pos" in msg
        assert "pip" in msg


# =============================================================================
# Gene annotation validation edge cases
# =============================================================================


class TestGenesValidationEdgeCases:
    """Edge cases for the gene annotation load-time contract."""

    def test_valid_genes_passes(self):
        """A minimal gene annotation frame satisfies the contract."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                "end": [2000],
                "gene_name": ["GENE1"],
            }
        )
        check(df, spec(Family.GENES, Tier.LOAD))

    def test_start_at_zero_passes(self):
        """Start position of 0 is valid (only negative is rejected)."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [0],
                "end": [1000],
                "gene_name": ["GENE1"],
            }
        )
        check(df, spec(Family.GENES, Tier.LOAD))

    def test_non_numeric_start_fails(self):
        """Non-numeric start column is rejected."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": ["abc"],
                "end": [1000],
                "gene_name": ["GENE1"],
            }
        )
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.GENES, Tier.LOAD))

    def test_non_numeric_end_fails(self):
        """Non-numeric end column is rejected."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                "end": ["xyz"],
                "gene_name": ["GENE1"],
            }
        )
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.GENES, Tier.LOAD))

    def test_equal_start_end_passes(self):
        """Gene with start == end is valid (e.g. single-base feature)."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                "end": [1000],
                "gene_name": ["POINT"],
            }
        )
        check(df, spec(Family.GENES, Tier.LOAD))

    def test_missing_all_required_columns(self):
        """All four missing columns reported."""
        df = pd.DataFrame({"extra": [1]})
        with pytest.raises(LoaderValidationError) as exc_info:
            check(df, spec(Family.GENES, Tier.LOAD))
        msg = str(exc_info.value)
        for col in ["chr", "start", "end", "gene_name"]:
            assert col in msg

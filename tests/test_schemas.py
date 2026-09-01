"""Tests for schemas.py validation functions.

Tests edge cases and error paths for GWAS, eQTL, fine-mapping, and gene
annotation DataFrame validation, and file path validation.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import LoaderValidationError
from pylocuszoom.schemas import (
    validate_eqtl_dataframe,
    validate_finemapping_dataframe,
    validate_genes_dataframe,
    validate_gwas_dataframe,
)

# =============================================================================
# GWAS validation edge cases
# =============================================================================


class TestGWASValidationEdgeCases:
    """Edge cases for validate_gwas_dataframe beyond basic happy/error paths."""

    def test_valid_gwas_returns_same_dataframe(self):
        """Return value is the same DataFrame object that was passed in."""
        df = pd.DataFrame({"ps": [100], "p_wald": [0.5], "rs": ["rs1"]})
        result = validate_gwas_dataframe(df)
        assert result is df

    def test_custom_column_names_validated(self):
        """Custom pos_col and p_col names are respected during validation."""
        df = pd.DataFrame({"position": [100], "pvalue": [0.5]})
        result = validate_gwas_dataframe(df, pos_col="position", p_col="pvalue")
        assert result is df

    def test_custom_column_names_missing_raises(self):
        """Missing custom column names raise LoaderValidationError."""
        df = pd.DataFrame({"ps": [100], "p_wald": [0.5]})
        with pytest.raises(LoaderValidationError, match=r"Missing columns.*'pos'"):
            validate_gwas_dataframe(df, pos_col="pos", p_col="p_wald")

    def test_both_columns_missing_reports_both(self):
        """Error message mentions both missing columns when both are absent."""
        df = pd.DataFrame({"chr": [1], "rs": ["rs1"]})
        with pytest.raises(LoaderValidationError, match="ps") as exc_info:
            validate_gwas_dataframe(df)
        assert "p_wald" in str(exc_info.value)

    def test_zero_position_fails(self):
        """Position of exactly zero is rejected (non-positive)."""
        df = pd.DataFrame({"ps": [0], "p_wald": [0.5]})
        with pytest.raises(LoaderValidationError, match=r"'ps': 1 values <= 0"):
            validate_gwas_dataframe(df)

    def test_pvalue_exactly_one_passes(self):
        """P-value of exactly 1.0 is valid (range is (0, 1])."""
        df = pd.DataFrame({"ps": [100], "p_wald": [1.0]})
        result = validate_gwas_dataframe(df)
        assert result is df

    def test_nan_pvalue_fails(self):
        """NaN p-values are caught."""
        df = pd.DataFrame({"ps": [100, 200], "p_wald": [0.05, float("nan")]})
        with pytest.raises(LoaderValidationError, match=r"'p_wald' has 1 null values"):
            validate_gwas_dataframe(df)

    def test_multiple_numeric_errors_reported(self):
        """Multiple validation errors accumulated and reported at once."""
        df = pd.DataFrame(
            {
                "ps": ["chr1:100"],
                "p_wald": ["not_a_number"],
            }  # non-numeric  # non-numeric
        )
        with pytest.raises(LoaderValidationError, match="must be numeric") as exc_info:
            validate_gwas_dataframe(df)
        msg = str(exc_info.value)
        assert "ps" in msg
        assert "p_wald" in msg

    def test_single_row_valid(self):
        """Single-row DataFrame passes validation."""
        df = pd.DataFrame({"ps": [42], "p_wald": [1e-10]})
        assert validate_gwas_dataframe(df) is df

    def test_extra_columns_ignored(self):
        """Extra columns beyond required ones don't cause failures."""
        df = pd.DataFrame(
            {"ps": [100], "p_wald": [0.5], "rs": ["rs1"], "beta": [0.3], "se": [0.1]}
        )
        assert validate_gwas_dataframe(df) is df


# =============================================================================
# eQTL validation edge cases
# =============================================================================


class TestEQTLValidationEdgeCases:
    """Edge cases for validate_eqtl_dataframe."""

    def test_valid_eqtl_returns_dataframe(self):
        """Valid eQTL DataFrame is returned unchanged."""
        df = pd.DataFrame({"pos": [100], "p_value": [0.05], "gene": ["BRCA1"]})
        assert validate_eqtl_dataframe(df) is df

    def test_non_numeric_pos_fails(self):
        """Non-numeric position column is rejected."""
        df = pd.DataFrame({"pos": ["chr1:100"], "p_value": [0.05], "gene": ["BRCA1"]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            validate_eqtl_dataframe(df)

    def test_non_numeric_pvalue_fails(self):
        """Non-numeric p_value column is rejected."""
        df = pd.DataFrame({"pos": [100], "p_value": ["significant"], "gene": ["BRCA1"]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            validate_eqtl_dataframe(df)

    def test_zero_pvalue_fails(self):
        """P-value of 0 is rejected (outside (0, 1])."""
        df = pd.DataFrame({"pos": [100], "p_value": [0.0], "gene": ["BRCA1"]})
        with pytest.raises(LoaderValidationError, match=r"'p_value': 1 values <= 0"):
            validate_eqtl_dataframe(df)

    def test_missing_multiple_columns(self):
        """All missing required columns listed in error."""
        df = pd.DataFrame({"extra": [1]})
        with pytest.raises(LoaderValidationError) as exc_info:
            validate_eqtl_dataframe(df)
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
            validate_eqtl_dataframe(df)


# =============================================================================
# Fine-mapping validation edge cases
# =============================================================================


class TestFinemappingValidationEdgeCases:
    """Edge cases for validate_finemapping_dataframe."""

    def test_pip_zero_and_one_valid(self):
        """PIP of 0.0 and 1.0 are both valid (range [0, 1])."""
        df = pd.DataFrame({"pos": [100, 200], "pip": [0.0, 1.0]})
        assert validate_finemapping_dataframe(df) is df

    def test_non_numeric_pip_fails(self):
        """Non-numeric pip column is rejected."""
        df = pd.DataFrame({"pos": [100], "pip": ["high"]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            validate_finemapping_dataframe(df)

    def test_non_numeric_pos_fails(self):
        """Non-numeric pos column is rejected."""
        df = pd.DataFrame({"pos": ["abc"], "pip": [0.5]})
        with pytest.raises(LoaderValidationError, match="must be numeric"):
            validate_finemapping_dataframe(df)

    def test_missing_both_columns(self):
        """Both pos and pip missing reported."""
        df = pd.DataFrame({"other": [1]})
        with pytest.raises(LoaderValidationError) as exc_info:
            validate_finemapping_dataframe(df)
        msg = str(exc_info.value)
        assert "pos" in msg
        assert "pip" in msg


# =============================================================================
# Gene annotation validation edge cases
# =============================================================================


class TestGenesValidationEdgeCases:
    """Edge cases for validate_genes_dataframe."""

    def test_valid_genes_returns_dataframe(self):
        """Valid genes DataFrame is returned unchanged."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                "end": [2000],
                "gene_name": ["GENE1"],
            }
        )
        assert validate_genes_dataframe(df) is df

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
        assert validate_genes_dataframe(df) is df

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
            validate_genes_dataframe(df)

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
            validate_genes_dataframe(df)

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
        assert validate_genes_dataframe(df) is df

    def test_missing_all_required_columns(self):
        """All four missing columns reported."""
        df = pd.DataFrame({"extra": [1]})
        with pytest.raises(LoaderValidationError) as exc_info:
            validate_genes_dataframe(df)
        msg = str(exc_info.value)
        for col in ["chr", "start", "end", "gene_name"]:
            assert col in msg

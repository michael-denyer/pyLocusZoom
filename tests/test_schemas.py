"""Tests for the load-time contracts in schemas.py.

One class per family. Each owns the required columns, the value domains and
the dtype rules that ``check(df, spec(family, Tier.LOAD))`` enforces.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import LoaderValidationError
from pylocuszoom.schemas import Family, Tier, spec
from pylocuszoom.validation import check


class TestGWASValidation:
    """The GWAS load-time contract: required columns, domains, and dtypes."""

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

    def test_valid_multi_row_gwas_passes(self):
        """A multi-row GWAS frame satisfies the contract."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000, 1002000],
                "p_wald": [0.01, 0.001, 1e-8],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )

        check(df, spec(Family.GWAS, Tier.LOAD))

    def test_missing_position_column_fails(self):
        """A frame without the position column is rejected."""
        df = pd.DataFrame(
            {
                "p_wald": [0.01, 0.001],
                "rs": ["rs1", "rs2"],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing columns"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_missing_pvalue_column_fails(self):
        """A frame without the p-value column is rejected."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "rs": ["rs1", "rs2"],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing columns"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_negative_position_fails(self):
        """Negative positions are rejected."""
        df = pd.DataFrame(
            {
                "ps": [-1000, 1001000],
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values <= 0"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_pvalue_out_of_range_fails(self):
        """P-values above 1 are rejected."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": [0.01, 1.5],  # 1.5 is out of range
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values > 1"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_zero_pvalue_fails(self):
        """A p-value of 0 is rejected."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": [0.0, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values <= 0"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_nan_position_fails(self):
        """NaN positions are rejected."""
        df = pd.DataFrame(
            {
                "ps": [1000000, None],
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"null values"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_non_numeric_pvalue_fails(self):
        """Non-numeric p-values are rejected."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": ["0.01", "significant"],  # Strings, not numbers
            }
        )

        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.GWAS, Tier.LOAD))

    def test_non_numeric_position_fails(self):
        """Non-numeric positions are rejected."""
        df = pd.DataFrame(
            {
                "ps": ["chr1:1000", "chr1:2000"],  # Strings, not numbers
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match="must be numeric"):
            check(df, spec(Family.GWAS, Tier.LOAD))


class TestEQTLValidation:
    """The eQTL load-time contract: required columns, domains, and dtypes."""

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

    def test_valid_eqtl_with_effect_passes(self):
        """An eQTL frame carrying signed effect sizes satisfies the contract."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
                "gene": ["BRCA1", "BRCA1"],
                "effect": [0.5, -0.3],
            }
        )

        check(df, spec(Family.EQTL, Tier.LOAD))

    def test_missing_gene_column_fails(self):
        """A frame without the gene column is rejected."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing columns"):
            check(df, spec(Family.EQTL, Tier.LOAD))


class TestFinemappingValidation:
    """The fine-mapping load-time contract: required columns and PIP domain."""

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

    def test_valid_finemapping_passes(self):
        """A fine-mapping frame with credible-set labels satisfies the contract."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000, 1002000],
                "pip": [0.85, 0.12, 0.03],
                "cs": [1, 1, 0],
            }
        )

        check(df, spec(Family.FINEMAPPING, Tier.LOAD))

    def test_pip_out_of_range_fails(self):
        """A PIP above 1 is rejected."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "pip": [0.85, 1.5],  # 1.5 is out of range
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values > 1"):
            check(df, spec(Family.FINEMAPPING, Tier.LOAD))

    def test_negative_pip_fails(self):
        """A negative PIP is rejected."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "pip": [-0.1, 0.5],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values < 0"):
            check(df, spec(Family.FINEMAPPING, Tier.LOAD))


class TestGenesValidation:
    """The gene annotation load-time contract: required columns and coordinates."""

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

    def test_valid_multi_row_genes_passes(self):
        """A multi-row gene annotation frame satisfies the contract."""
        df = pd.DataFrame(
            {
                "chr": ["1", "1", "1"],
                "start": [1000000, 1050000, 1100000],
                "end": [1020000, 1080000, 1150000],
                "gene_name": ["GENE1", "GENE2", "GENE3"],
            }
        )

        check(df, spec(Family.GENES, Tier.LOAD))

    def test_end_before_start_fails(self):
        """A gene whose end precedes its start is rejected."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1020000],  # Start after end
                "end": [1000000],
                "gene_name": ["GENE1"],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"rows have start > end"):
            check(df, spec(Family.GENES, Tier.LOAD))

    def test_negative_start_fails(self):
        """A negative start coordinate is rejected."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [-1000],
                "end": [1000000],
                "gene_name": ["GENE1"],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values < 0"):
            check(df, spec(Family.GENES, Tier.LOAD))

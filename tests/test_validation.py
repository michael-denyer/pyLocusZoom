"""Tests for the ColumnSpec validation engine."""

import numpy as np
import pandas as pd
import pytest
from hypothesis import given

from pylocuszoom.exceptions import ValidationError
from pylocuszoom.schemas import validate_genes_df, validate_gwas_df
from pylocuszoom.validation import ColumnSpec, RangeRule, check
from tests.strategies import gwas_dataframes, pvalues, pvalues_invalid


class TestRequiredColumns:
    """Test the required-columns rule."""

    def test_all_columns_present(self):
        """No error when all required columns exist."""
        df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
        check(df, ColumnSpec(name="test_df", required=("a", "b")))

    def test_missing_single_column(self):
        """Error message includes missing and available columns."""
        df = pd.DataFrame({"a": [1]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", required=("a", "b")))

        error_msg = str(exc_info.value)
        assert "test_df validation failed" in error_msg
        assert "Missing columns: ['b']" in error_msg
        assert "Available: ['a']" in error_msg

    def test_missing_multiple_columns(self):
        """Error lists all missing columns."""
        df = pd.DataFrame({"a": [1]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", required=("a", "b", "c")))

        assert "Missing columns: ['b', 'c']" in str(exc_info.value)

    def test_empty_columns_list(self):
        """No error when no columns required."""
        df = pd.DataFrame({"a": [1]})
        check(df, ColumnSpec(name="test_df"))


class TestNumericColumns:
    """Test the numeric-dtype rule."""

    def test_numeric_int_column(self):
        """No error for integer column."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        check(df, ColumnSpec(name="test_df", numeric=("a",)))

    def test_numeric_float_column(self):
        """No error for float column."""
        df = pd.DataFrame({"a": [1.5, 2.5, 3.5]})
        check(df, ColumnSpec(name="test_df", numeric=("a",)))

    def test_non_numeric_string_column(self):
        """Error for string column."""
        df = pd.DataFrame({"a": ["x", "y", "z"]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", numeric=("a",)))

        error_msg = str(exc_info.value)
        assert "Column 'a' must be numeric" in error_msg
        assert "str" in error_msg or "object" in error_msg

    def test_multiple_non_numeric_columns(self):
        """Error lists all non-numeric columns."""
        df = pd.DataFrame({"a": ["x", "y"], "b": ["p", "q"], "c": [1, 2]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", numeric=("a", "b")))

        error_msg = str(exc_info.value)
        assert "Column 'a' must be numeric" in error_msg
        assert "Column 'b' must be numeric" in error_msg

    def test_skip_missing_columns(self):
        """Don't check dtype for columns that don't exist."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        check(df, ColumnSpec(name="test_df", numeric=("b",)))

    def test_non_numeric_column_short_circuits_its_range_rule(self):
        """A non-numeric column reports a dtype error, never a TypeError."""
        df = pd.DataFrame({"a": ["x", "y", "z"]})

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test_df",
                    numeric=("a",),
                    ranges=(RangeRule("a", min_val=0, max_val=1),),
                ),
            )

        error_msg = str(exc_info.value)
        assert "Column 'a' must be numeric" in error_msg
        assert "values" not in error_msg


class TestRangeRules:
    """Test the numeric-range rule."""

    def test_values_within_range(self):
        """No error when all values within bounds."""
        df = pd.DataFrame({"p": [0.1, 0.5, 0.9]})
        check(
            df,
            ColumnSpec(name="test_df", ranges=(RangeRule("p", min_val=0, max_val=1),)),
        )

    def test_values_exceed_max(self):
        """Error when values exceed max."""
        df = pd.DataFrame({"p": [0.5, 1.5, 2.0]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", ranges=(RangeRule("p", max_val=1),)))

        assert "Column 'p': 2 values > 1" in str(exc_info.value)

    def test_values_below_min(self):
        """Error when values below min."""
        df = pd.DataFrame({"p": [-1.0, 0.0, 0.5]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", ranges=(RangeRule("p", min_val=0),)))

        assert "Column 'p': 1 values < 0" in str(exc_info.value)

    def test_nan_values_ignored_in_range_check(self):
        """NaN values should be skipped during range validation."""
        df = pd.DataFrame({"p": [0.1, np.nan, 0.5, np.nan, 0.9]})
        check(
            df,
            ColumnSpec(name="test_df", ranges=(RangeRule("p", min_val=0, max_val=1),)),
        )

    def test_all_nan_column_passes_range_check(self):
        """Column with all NaN values should pass range check."""
        df = pd.DataFrame({"p": [np.nan, np.nan, np.nan]})
        check(
            df,
            ColumnSpec(name="test_df", ranges=(RangeRule("p", min_val=0, max_val=1),)),
        )

    def test_exclusive_min(self):
        """Error when values equal to exclusive min."""
        df = pd.DataFrame({"p": [0.0, 0.5, 1.0]})

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test_df",
                    ranges=(RangeRule("p", min_val=0, exclusive_min=True),),
                ),
            )

        assert "Column 'p': 1 values <= 0" in str(exc_info.value)

    def test_exclusive_max(self):
        """Error when values equal to exclusive max."""
        df = pd.DataFrame({"p": [0.0, 0.5, 1.0]})

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test_df",
                    ranges=(RangeRule("p", max_val=1, exclusive_max=True),),
                ),
            )

        assert "Column 'p': 1 values >= 1" in str(exc_info.value)

    def test_min_and_max(self):
        """Check both bounds."""
        df = pd.DataFrame({"p": [-0.5, 0.5, 1.5]})

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test_df", ranges=(RangeRule("p", min_val=0, max_val=1),)
                ),
            )

        error_msg = str(exc_info.value)
        assert "1 values < 0" in error_msg
        assert "1 values > 1" in error_msg

    def test_skip_missing_column(self):
        """Don't check range for missing column."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        check(
            df,
            ColumnSpec(name="test_df", ranges=(RangeRule("b", min_val=0, max_val=10),)),
        )


class TestNotNullColumns:
    """Test the not-null rule."""

    def test_no_null_values(self):
        """No error when no nulls."""
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        check(df, ColumnSpec(name="test_df", not_null=("a", "b")))

    def test_nan_values(self):
        """Error when NaN present."""
        df = pd.DataFrame({"a": [1, np.nan, 3]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", not_null=("a",)))

        assert "Column 'a' has 1 null values" in str(exc_info.value)

    def test_multiple_null_values(self):
        """Report count of nulls."""
        df = pd.DataFrame({"a": [1, np.nan, np.nan, 4]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", not_null=("a",)))

        assert "Column 'a' has 2 null values" in str(exc_info.value)

    def test_none_values(self):
        """Error when None present."""
        df = pd.DataFrame({"a": [1, None, 3]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", not_null=("a",)))

        assert "Column 'a' has 1 null values" in str(exc_info.value)

    def test_multiple_columns_with_nulls(self):
        """Report nulls in multiple columns."""
        df = pd.DataFrame({"a": [1, np.nan], "b": [None, 2]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", not_null=("a", "b")))

        error_msg = str(exc_info.value)
        assert "Column 'a' has 1 null values" in error_msg
        assert "Column 'b' has 1 null values" in error_msg

    def test_skip_missing_column(self):
        """Don't check nulls for missing column."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        check(df, ColumnSpec(name="test_df", not_null=("b",)))


class TestNonEmpty:
    """Test the non-empty rule."""

    def test_empty_frame_rejected(self):
        """A frame with no rows is rejected before any column rule runs."""
        df = pd.DataFrame({"a": pd.Series(dtype="int64")})

        with pytest.raises(ValidationError, match="is empty"):
            check(df, ColumnSpec(name="test_df", non_empty=True))

    def test_empty_frame_allowed_when_rule_is_off(self):
        """Without non_empty, a row-less frame with the right columns passes."""
        df = pd.DataFrame({"a": pd.Series(dtype="int64")})
        check(df, ColumnSpec(name="test_df", required=("a",)))


class TestOrdering:
    """Test the row-wise ordering rule."""

    def test_inverted_rows_reported_with_count(self):
        """Rows where lower exceeds upper are counted in the error."""
        df = pd.DataFrame({"start": [10, 30, 50], "end": [20, 25, 40]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test_df", ordering=(("start", "end"),)))

        assert "2 rows have start > end" in str(exc_info.value)

    def test_equal_bounds_pass(self):
        """Equal lower and upper bounds are ordered, not inverted."""
        df = pd.DataFrame({"start": [10], "end": [10]})
        check(df, ColumnSpec(name="test_df", ordering=(("start", "end"),)))

    def test_skip_missing_column(self):
        """Ordering is not checked when a bound column is absent."""
        df = pd.DataFrame({"start": [10]})
        check(df, ColumnSpec(name="test_df", ordering=(("start", "end"),)))

    def test_skip_non_numeric_column(self):
        """A non-numeric bound reports its dtype, not an ordering comparison."""
        df = pd.DataFrame({"start": ["a"], "end": [10]})

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test_df",
                    numeric=("start",),
                    ordering=(("start", "end"),),
                ),
            )

        error_msg = str(exc_info.value)
        assert "must be numeric" in error_msg
        assert "rows have" not in error_msg


class TestPValueDomain:
    """Test the canonical (0, 1] p-value rule."""

    def test_pvalue_of_one_passes(self):
        """One is the closed upper bound of the p-value domain."""
        df = pd.DataFrame({"p": [1.0, 0.5]})
        check(df, ColumnSpec(name="test_df", pvalue="p"))

    def test_pvalue_of_zero_rejected(self):
        """Zero is excluded by the open lower bound."""
        df = pd.DataFrame({"p": [0.0, 0.5]})

        with pytest.raises(ValidationError, match=r"'p': 1 values <= 0"):
            check(df, ColumnSpec(name="test_df", pvalue="p"))

    def test_pvalue_above_one_rejected(self):
        """Anything above one is outside the domain."""
        df = pd.DataFrame({"p": [1.5, 0.5]})

        with pytest.raises(ValidationError, match=r"'p': 1 values > 1"):
            check(df, ColumnSpec(name="test_df", pvalue="p"))


class TestErrorAccumulation:
    """Test that multiple errors are accumulated and reported together."""

    def test_accumulate_multiple_errors(self):
        """All errors should be reported in single ValidationError."""
        df = pd.DataFrame(
            {
                "a": [1, 2, 3],
                "b": ["x", "y", "z"],
                "c": [0.5, 1.5, 2.5],
            }
        )

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test_df",
                    required=("a", "b", "c", "d"),
                    numeric=("b",),
                    ranges=(RangeRule("c", min_val=0, max_val=1),),
                ),
            )

        error_msg = str(exc_info.value)
        assert "Missing columns: ['d']" in error_msg
        assert "Column 'b' must be numeric" in error_msg
        assert "Column 'c': 2 values > 1" in error_msg

    def test_no_errors_accumulated(self):
        """check() succeeds when no errors."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        check(
            df,
            ColumnSpec(
                name="test_df",
                required=("a",),
                numeric=("a",),
                not_null=("a",),
                ranges=(RangeRule("a", min_val=0, max_val=10),),
            ),
        )


class TestCustomName:
    """Test that the spec name appears in error messages."""

    def test_custom_name_in_error(self):
        """Error message should include the spec name."""
        df = pd.DataFrame({"a": [1]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="gwas_df", required=("b",)))

        assert "gwas_df validation failed" in str(exc_info.value)


class TestErrorClass:
    """Test that the spec's concrete exception type is raised."""

    def test_spec_error_class_is_raised(self):
        """The raised type is spec.error_class, not the base ValidationError."""
        from pylocuszoom.exceptions import LoaderValidationError

        df = pd.DataFrame({"a": [1]})

        with pytest.raises(LoaderValidationError):
            check(
                df,
                ColumnSpec(
                    name="test_df",
                    required=("b",),
                    error_class=LoaderValidationError,
                ),
            )


class TestErrorPathCoverage:
    """Test error message content for validation failures.

    These tests ensure validation errors provide actionable information
    for debugging. Error messages should include context needed to fix
    the issue without additional investigation.
    """

    def test_missing_columns_error_shows_available(self):
        """Error message should list both missing and available columns."""
        df = pd.DataFrame({"a": [1], "b": [2]})

        with pytest.raises(ValidationError, match=r"Available.*\['a', 'b'\]"):
            check(df, ColumnSpec(name="test", required=("x", "y")))

    def test_non_numeric_column_shows_dtype(self):
        """Non-numeric column error should show actual dtype."""
        df = pd.DataFrame({"val": ["a", "b", "c"]})

        with pytest.raises(ValidationError, match=r"must be numeric, got (object|str)"):
            check(df, ColumnSpec(name="test", numeric=("val",)))

    def test_range_out_of_bounds_shows_count(self):
        """Out-of-range error should report count and bound violated."""
        df = pd.DataFrame({"pval": [0.1, -0.5, 1.5]})

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test", ranges=(RangeRule("pval", min_val=0, max_val=1),)
                ),
            )

        error_msg = str(exc_info.value)
        assert "1 values < 0" in error_msg
        assert "1 values > 1" in error_msg

    def test_error_includes_dataframe_name(self):
        """Error header should include DataFrame name for context."""
        df = pd.DataFrame({"x": [1]})

        with pytest.raises(ValidationError, match=r"GWAS DataFrame validation failed"):
            check(df, ColumnSpec(name="GWAS DataFrame", required=("pos",)))

    def test_multiple_errors_all_reported(self):
        """All validation errors should be accumulated and reported."""
        df = pd.DataFrame(
            {
                "numeric_col": [1, 2, 3],
                "string_col": ["a", "b", "c"],
                "range_col": [0.5, 1.5, 2.5],
            }
        )

        with pytest.raises(ValidationError) as exc_info:
            check(
                df,
                ColumnSpec(
                    name="test",
                    required=("missing_col",),
                    numeric=("string_col",),
                    ranges=(RangeRule("range_col", max_val=1),),
                ),
            )

        error_msg = str(exc_info.value)
        assert "Missing columns: ['missing_col']" in error_msg
        assert "Column 'string_col' must be numeric" in error_msg
        assert "2 values > 1" in error_msg

    def test_null_check_reports_exact_count(self):
        """Null check should report exact count of nulls."""
        df = pd.DataFrame({"val": [1, np.nan, np.nan, np.nan, 5]})

        with pytest.raises(ValidationError, match=r"has 3 null values"):
            check(df, ColumnSpec(name="test", not_null=("val",)))

    def test_error_message_formatted_as_list(self):
        """Error message should format issues as readable list."""
        df = pd.DataFrame({"a": ["x"]})

        with pytest.raises(ValidationError) as exc_info:
            check(df, ColumnSpec(name="test", required=("b",), numeric=("a",)))

        error_msg = str(exc_info.value)
        assert "  - Missing columns:" in error_msg
        assert "  - Column 'a' must be numeric" in error_msg


class TestEQTLValidation:
    """Tests for eQTL-specific validation.

    Bug fix: pyLocusZoom-7a5
    validate_eqtl_df should enforce numeric p_value column.
    """

    def test_validate_eqtl_df_requires_numeric_pvalue(self):
        """eQTL validation should fail for non-numeric p_value."""
        from pylocuszoom.eqtl import validate_eqtl_df
        from pylocuszoom.exceptions import EQTLValidationError

        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "p_value": ["0.01", "0.05", "0.001"],  # Strings, not floats
            }
        )

        with pytest.raises(EQTLValidationError, match="numeric"):
            validate_eqtl_df(df, pos_col="pos", p_col="p_value")

    def test_validate_eqtl_df_accepts_numeric_pvalue(self):
        """eQTL validation should pass for numeric p_value."""
        from pylocuszoom.eqtl import validate_eqtl_df

        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "p_value": [0.01, 0.05, 0.001],  # Numeric
            }
        )

        # Should not raise
        validate_eqtl_df(df, pos_col="pos", p_col="p_value")


# =============================================================================
# Property-Based Tests (Hypothesis)
# =============================================================================


class TestValidationProperties:
    """Property-based tests for DataFrame validation."""

    @given(gwas_dataframes(min_snps=1, max_snps=50))
    def test_valid_gwas_df_passes_validation(self, df):
        """Valid GWAS DataFrames should always pass validation."""
        check(
            df,
            ColumnSpec(
                name="gwas_df",
                required=("rs", "chr", "ps", "p_wald"),
                numeric=("ps", "p_wald"),
            ),
        )

    @given(gwas_dataframes(min_snps=1, max_snps=50))
    def test_dropping_required_column_fails(self, df):
        """Dropping any required column should fail validation."""
        required = ("rs", "chr", "ps", "p_wald")
        for col in required:
            invalid_df = df.drop(columns=[col])
            with pytest.raises(ValidationError):
                check(invalid_df, ColumnSpec(name="test", required=required))


class TestRangeValidationProperties:
    """Property tests for range validation."""

    @given(pvalues(allow_zero=False, allow_one=True))
    def test_valid_pvalue_in_range(self, p):
        """Valid p-values should pass range check for (0, 1]."""
        df = pd.DataFrame({"p": [p]})
        check(df, ColumnSpec(name="test", pvalue="p"))

    @given(pvalues_invalid().filter(lambda p: not np.isnan(p)))
    def test_invalid_pvalue_detected(self, p):
        """Invalid p-values should fail range validation."""
        df = pd.DataFrame({"p": [p]})
        with pytest.raises(ValidationError):
            check(df, ColumnSpec(name="test", pvalue="p"))


class TestValidateGwasDf:
    """Tests for validate_gwas_df function."""

    def test_valid_gwas_passes(self):
        """Valid GWAS DataFrame passes."""
        df = pd.DataFrame({"ps": [1000], "p_wald": [0.01]})
        validate_gwas_df(df)  # Should not raise

    def test_custom_column_names(self):
        """Custom column names work."""
        df = pd.DataFrame({"pos": [1000], "pval": [0.01]})
        validate_gwas_df(df, pos_col="pos", p_col="pval")  # Should not raise

    def test_missing_position_raises(self):
        """Missing position column raises error."""
        df = pd.DataFrame({"p_wald": [0.01]})

        with pytest.raises(ValidationError):
            validate_gwas_df(df)

    def test_with_rs_col(self):
        """Including rs_col validates that column too."""
        df = pd.DataFrame({"ps": [1000], "p_wald": [0.01], "rs": ["rs123"]})
        validate_gwas_df(df, rs_col="rs")  # Should not raise

    def test_missing_rs_col_when_required(self):
        """Missing rs_col when specified raises error."""
        df = pd.DataFrame({"ps": [1000], "p_wald": [0.01]})

        with pytest.raises(ValidationError):
            validate_gwas_df(df, rs_col="rs")


class TestValidateGenesDf:
    """Tests for validate_genes_df function."""

    def test_valid_genes_passes(self):
        """Valid genes DataFrame passes."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                "end": [2000],
                "gene_name": ["BRCA1"],
            }
        )
        validate_genes_df(df)  # Should not raise

    def test_missing_column_raises(self):
        """Missing required column raises error."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                # Missing "end" and "gene_name"
            }
        )

        with pytest.raises(ValidationError):
            validate_genes_df(df)

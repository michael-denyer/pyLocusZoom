"""Tests for DataFrameValidator builder class."""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.utils import ValidationError
from pylocuszoom.validation import DataFrameValidator


class TestRequireColumns:
    """Test require_columns() method."""

    def test_all_columns_present(self):
        """No error when all required columns exist."""
        df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_columns(["a", "b"]).validate()
        # Should not raise

    def test_missing_single_column(self):
        """Error message includes missing and available columns."""
        df = pd.DataFrame({"a": [1]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_columns(["a", "b"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "test_df validation failed" in error_msg
        assert "Missing columns: ['b']" in error_msg
        assert "Available: ['a']" in error_msg

    def test_missing_multiple_columns(self):
        """Error lists all missing columns."""
        df = pd.DataFrame({"a": [1]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_columns(["a", "b", "c"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Missing columns: ['b', 'c']" in error_msg

    def test_empty_columns_list(self):
        """No error when no columns required."""
        df = pd.DataFrame({"a": [1]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_columns([]).validate()
        # Should not raise


class TestRequireNumeric:
    """Test require_numeric() method."""

    def test_numeric_int_column(self):
        """No error for integer column."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_numeric(["a"]).validate()
        # Should not raise

    def test_numeric_float_column(self):
        """No error for float column."""
        df = pd.DataFrame({"a": [1.5, 2.5, 3.5]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_numeric(["a"]).validate()
        # Should not raise

    def test_non_numeric_string_column(self):
        """Error for string column."""
        df = pd.DataFrame({"a": ["x", "y", "z"]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_numeric(["a"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'a' must be numeric" in error_msg
        assert "str" in error_msg or "object" in error_msg

    def test_multiple_non_numeric_columns(self):
        """Error lists all non-numeric columns."""
        df = pd.DataFrame({"a": ["x", "y"], "b": ["p", "q"], "c": [1, 2]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_numeric(["a", "b"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'a' must be numeric" in error_msg
        assert "Column 'b' must be numeric" in error_msg

    def test_skip_missing_columns(self):
        """Don't check dtype for columns that don't exist."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        validator = DataFrameValidator(df, name="test_df")
        # Only require_numeric should complain about missing column
        validator.require_numeric(["b"])
        # This should NOT raise about dtype, only about missing column
        # if require_columns was called
        validator.validate()
        # Should not raise - missing columns handled by require_columns


class TestRequireRange:
    """Test require_range() method."""

    def test_values_within_range(self):
        """No error when all values within bounds."""
        df = pd.DataFrame({"p": [0.1, 0.5, 0.9]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("p", min_val=0, max_val=1).validate()
        # Should not raise

    def test_values_exceed_max(self):
        """Error when values exceed max."""
        df = pd.DataFrame({"p": [0.5, 1.5, 2.0]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("p", max_val=1)

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'p': 2 values > 1" in error_msg

    def test_values_below_min(self):
        """Error when values below min."""
        df = pd.DataFrame({"p": [-1.0, 0.0, 0.5]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("p", min_val=0)

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'p': 1 values < 0" in error_msg

    def test_exclusive_min(self):
        """Error when values equal to exclusive min."""
        df = pd.DataFrame({"p": [0.0, 0.5, 1.0]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("p", min_val=0, exclusive_min=True)

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'p': 1 values <= 0" in error_msg

    def test_exclusive_max(self):
        """Error when values equal to exclusive max."""
        df = pd.DataFrame({"p": [0.0, 0.5, 1.0]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("p", max_val=1, exclusive_max=True)

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'p': 1 values >= 1" in error_msg

    def test_min_and_max(self):
        """Check both bounds."""
        df = pd.DataFrame({"p": [-0.5, 0.5, 1.5]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("p", min_val=0, max_val=1)

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        # Should report both violations
        assert "1 values < 0" in error_msg
        assert "1 values > 1" in error_msg

    def test_skip_missing_column(self):
        """Don't check range for missing column."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_range("b", min_val=0, max_val=10)
        validator.validate()
        # Should not raise - missing columns handled separately


class TestRequireNotNull:
    """Test require_not_null() method."""

    def test_no_null_values(self):
        """No error when no nulls."""
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_not_null(["a", "b"]).validate()
        # Should not raise

    def test_nan_values(self):
        """Error when NaN present."""
        df = pd.DataFrame({"a": [1, np.nan, 3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_not_null(["a"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'a' has 1 null values" in error_msg

    def test_multiple_null_values(self):
        """Report count of nulls."""
        df = pd.DataFrame({"a": [1, np.nan, np.nan, 4]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_not_null(["a"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'a' has 2 null values" in error_msg

    def test_none_values(self):
        """Error when None present."""
        df = pd.DataFrame({"a": [1, None, 3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_not_null(["a"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'a' has 1 null values" in error_msg

    def test_multiple_columns_with_nulls(self):
        """Report nulls in multiple columns."""
        df = pd.DataFrame({"a": [1, np.nan], "b": [None, 2]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_not_null(["a", "b"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        assert "Column 'a' has 1 null values" in error_msg
        assert "Column 'b' has 1 null values" in error_msg

    def test_skip_missing_column(self):
        """Don't check nulls for missing column."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_not_null(["b"])
        validator.validate()
        # Should not raise - missing columns handled separately


class TestMethodChaining:
    """Test that methods return self for chaining."""

    def test_chaining(self):
        """All methods should return self."""
        df = pd.DataFrame({"a": [1, 2, 3], "b": [0.1, 0.5, 0.9]})
        validator = DataFrameValidator(df, name="test_df")

        # Chain all methods
        result = (
            validator.require_columns(["a", "b"])
            .require_numeric(["a", "b"])
            .require_range("a", min_val=0, max_val=10)
            .require_not_null(["a", "b"])
        )

        # Should return the validator instance
        assert result is validator

        # validate() should return None
        assert result.validate() is None


class TestErrorAccumulation:
    """Test that multiple errors are accumulated and reported together."""

    def test_accumulate_multiple_errors(self):
        """All errors should be reported in single ValidationError."""
        df = pd.DataFrame(
            {
                "a": [1, 2, 3],
                "b": ["x", "y", "z"],  # Not numeric
                "c": [0.5, 1.5, 2.5],  # Out of range
            }
        )
        validator = DataFrameValidator(df, name="test_df")
        validator.require_columns(["a", "b", "c", "d"])  # Missing 'd'
        validator.require_numeric(["b"])  # Wrong type
        validator.require_range("c", min_val=0, max_val=1)  # Out of range

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        error_msg = str(exc_info.value)
        # All three errors should be present
        assert "Missing columns: ['d']" in error_msg
        assert "Column 'b' must be numeric" in error_msg
        assert "Column 'c': 2 values > 1" in error_msg

    def test_no_errors_accumulated(self):
        """validate() succeeds when no errors."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        validator = DataFrameValidator(df, name="test_df")
        validator.require_columns(["a"])
        validator.require_numeric(["a"])
        validator.require_range("a", min_val=0, max_val=10)
        validator.require_not_null(["a"])
        validator.validate()
        # Should not raise


class TestCustomName:
    """Test that custom name appears in error messages."""

    def test_custom_name_in_error(self):
        """Error message should include custom name."""
        df = pd.DataFrame({"a": [1]})
        validator = DataFrameValidator(df, name="gwas_df")
        validator.require_columns(["b"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        assert "gwas_df validation failed" in str(exc_info.value)

    def test_default_name(self):
        """Default name should be 'DataFrame'."""
        df = pd.DataFrame({"a": [1]})
        validator = DataFrameValidator(df)
        validator.require_columns(["b"])

        with pytest.raises(ValidationError) as exc_info:
            validator.validate()

        assert "DataFrame validation failed" in str(exc_info.value)


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
            DataFrameValidator(df, "test").require_columns(["x", "y"]).validate()

    def test_require_numeric_non_numeric_column_shows_dtype(self):
        """Non-numeric column error should show actual dtype."""
        df = pd.DataFrame({"val": ["a", "b", "c"]})

        # dtype may be 'object' or 'str' depending on pandas version
        with pytest.raises(ValidationError, match=r"must be numeric, got (object|str)"):
            DataFrameValidator(df, "test").require_numeric(["val"]).validate()

    def test_require_range_out_of_bounds_shows_count(self):
        """Out-of-range error should report count and bound violated."""
        df = pd.DataFrame({"pval": [0.1, -0.5, 1.5]})

        with pytest.raises(ValidationError) as exc_info:
            DataFrameValidator(df, "test").require_range("pval", 0, 1).validate()

        error_msg = str(exc_info.value)
        # Should report both bound violations with counts
        assert "1 values < 0" in error_msg
        assert "1 values > 1" in error_msg

    def test_error_includes_dataframe_name(self):
        """Error header should include DataFrame name for context."""
        df = pd.DataFrame({"x": [1]})

        with pytest.raises(ValidationError, match=r"GWAS DataFrame validation failed"):
            DataFrameValidator(df, "GWAS DataFrame").require_columns(["pos"]).validate()

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
            (
                DataFrameValidator(df, "test")
                .require_columns(["missing_col"])
                .require_numeric(["string_col"])
                .require_range("range_col", max_val=1)
                .validate()
            )

        error_msg = str(exc_info.value)
        # All three issues reported
        assert "Missing columns: ['missing_col']" in error_msg
        assert "Column 'string_col' must be numeric" in error_msg
        assert "2 values > 1" in error_msg

    def test_null_check_reports_exact_count(self):
        """Null check should report exact count of nulls."""
        df = pd.DataFrame({"val": [1, np.nan, np.nan, np.nan, 5]})

        with pytest.raises(ValidationError, match=r"has 3 null values"):
            DataFrameValidator(df, "test").require_not_null(["val"]).validate()

    def test_error_message_formatted_as_list(self):
        """Error message should format issues as readable list."""
        df = pd.DataFrame({"a": ["x"]})

        with pytest.raises(ValidationError) as exc_info:
            (
                DataFrameValidator(df, "test")
                .require_columns(["b"])
                .require_numeric(["a"])
                .validate()
            )

        error_msg = str(exc_info.value)
        # Each error on its own line with bullet
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

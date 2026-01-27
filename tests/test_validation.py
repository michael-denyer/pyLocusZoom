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

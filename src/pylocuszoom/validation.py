"""DataFrame validation builder for pyLocusZoom.

Provides a fluent API for validating pandas DataFrames with composable
validation rules. Accumulates all validation errors before raising.
"""

from typing import List, Optional

import pandas as pd
from pandas.api.types import is_numeric_dtype

from .utils import ValidationError


class DataFrameValidator:
    """Builder for composable DataFrame validation.

    Validates DataFrames with method chaining and accumulates all errors
    before raising. This enables clear, readable validation code with
    comprehensive error messages.

    Example:
        >>> validator = DataFrameValidator(df, name="gwas_df")
        >>> validator.require_columns(["chr", "pos", "p"])
        ...     .require_numeric(["pos", "p"])
        ...     .require_range("p", min_val=0, max_val=1)
        ...     .validate()
    """

    def __init__(self, df: pd.DataFrame, name: str = "DataFrame"):
        """Initialize validator.

        Args:
            df: DataFrame to validate.
            name: Name for error messages (e.g., "gwas_df", "genes_df").
        """
        self._df = df
        self._name = name
        self._errors: List[str] = []

    def require_columns(self, columns: List[str]) -> "DataFrameValidator":
        """Check that required columns exist in DataFrame.

        Args:
            columns: List of required column names.

        Returns:
            Self for method chaining.
        """
        if not columns:
            return self

        missing = [col for col in columns if col not in self._df.columns]
        if missing:
            available = list(self._df.columns)
            self._errors.append(f"Missing columns: {missing}. Available: {available}")

        return self

    def require_numeric(self, columns: List[str]) -> "DataFrameValidator":
        """Check that columns have numeric dtype.

        Skips columns that don't exist (checked separately by require_columns).

        Args:
            columns: List of column names that should be numeric.

        Returns:
            Self for method chaining.
        """
        for col in columns:
            # Skip missing columns - let require_columns handle that
            if col not in self._df.columns:
                continue

            if not is_numeric_dtype(self._df[col]):
                actual_dtype = self._df[col].dtype
                self._errors.append(
                    f"Column '{col}' must be numeric, got {actual_dtype}"
                )

        return self

    def require_range(
        self,
        column: str,
        min_val: Optional[float] = None,
        max_val: Optional[float] = None,
        exclusive_min: bool = False,
        exclusive_max: bool = False,
    ) -> "DataFrameValidator":
        """Check that column values are within specified range.

        Args:
            column: Column name to check.
            min_val: Minimum allowed value (inclusive by default).
            max_val: Maximum allowed value (inclusive by default).
            exclusive_min: If True, minimum is exclusive (values must be > min_val).
            exclusive_max: If True, maximum is exclusive (values must be < max_val).

        Returns:
            Self for method chaining.
        """
        # Skip missing columns
        if column not in self._df.columns:
            return self

        col_data = self._df[column]

        # Check minimum bound
        if min_val is not None:
            if exclusive_min:
                invalid_count = (col_data <= min_val).sum()
                if invalid_count > 0:
                    self._errors.append(
                        f"Column '{column}': {invalid_count} values <= {min_val}"
                    )
            else:
                invalid_count = (col_data < min_val).sum()
                if invalid_count > 0:
                    self._errors.append(
                        f"Column '{column}': {invalid_count} values < {min_val}"
                    )

        # Check maximum bound
        if max_val is not None:
            if exclusive_max:
                invalid_count = (col_data >= max_val).sum()
                if invalid_count > 0:
                    self._errors.append(
                        f"Column '{column}': {invalid_count} values >= {max_val}"
                    )
            else:
                invalid_count = (col_data > max_val).sum()
                if invalid_count > 0:
                    self._errors.append(
                        f"Column '{column}': {invalid_count} values > {max_val}"
                    )

        return self

    def require_not_null(self, columns: List[str]) -> "DataFrameValidator":
        """Check that columns have no null (NaN or None) values.

        Args:
            columns: List of column names to check for nulls.

        Returns:
            Self for method chaining.
        """
        for col in columns:
            # Skip missing columns
            if col not in self._df.columns:
                continue

            null_count = self._df[col].isna().sum()
            if null_count > 0:
                self._errors.append(f"Column '{col}' has {null_count} null values")

        return self

    def validate(self) -> None:
        """Raise ValidationError if any validation rules failed.

        Raises:
            ValidationError: If any validation errors were accumulated.
                Error message includes all accumulated errors.
        """
        if self._errors:
            error_msg = f"{self._name} validation failed:\n"
            error_msg += "\n".join(f"  - {error}" for error in self._errors)
            raise ValidationError(error_msg)

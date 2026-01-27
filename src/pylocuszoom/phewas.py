"""PheWAS data validation and preparation.

Validates and prepares phenome-wide association study data for plotting.
"""

import pandas as pd

from .utils import ValidationError


def validate_phewas_df(
    df: pd.DataFrame,
    phenotype_col: str = "phenotype",
    p_col: str = "p_value",
    category_col: str = "category",
) -> None:
    """Validate PheWAS DataFrame has required columns and types.

    Args:
        df: PheWAS results DataFrame.
        phenotype_col: Column name for phenotype names.
        p_col: Column name for p-values.
        category_col: Column name for phenotype categories (optional).

    Raises:
        ValidationError: If required columns are missing or have invalid types.
    """
    errors = []

    # Check required columns exist
    required = [phenotype_col, p_col]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValidationError(
            f"PheWAS DataFrame missing required columns: {missing}. "
            f"Required: {required}. Found: {list(df.columns)}"
        )

    # Check p-value column is numeric
    if not pd.api.types.is_numeric_dtype(df[p_col]):
        errors.append(f"Column '{p_col}' must be numeric, got {df[p_col].dtype}")

    # Check p-value range (0, 1]
    if pd.api.types.is_numeric_dtype(df[p_col]):
        invalid_p = ((df[p_col] <= 0) | (df[p_col] > 1)).sum()
        if invalid_p > 0:
            errors.append(
                f"Column '{p_col}' has {invalid_p} values outside range (0, 1]"
            )

    if errors:
        raise ValidationError("PheWAS validation failed:\n  - " + "\n  - ".join(errors))

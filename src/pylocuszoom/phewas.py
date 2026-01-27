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
    """Validate PheWAS DataFrame has required columns.

    Args:
        df: PheWAS results DataFrame.
        phenotype_col: Column name for phenotype names.
        p_col: Column name for p-values.
        category_col: Column name for phenotype categories (optional).

    Raises:
        ValidationError: If required columns are missing.
    """
    required = [phenotype_col, p_col]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValidationError(
            f"PheWAS DataFrame missing required columns: {missing}. "
            f"Required: {required}. Found: {list(df.columns)}"
        )

"""PheWAS data validation and preparation.

Validates and prepares phenome-wide association study data for plotting.
"""

import pandas as pd

from .validation import DataFrameValidator


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
    (
        DataFrameValidator(df, "PheWAS DataFrame")
        .require_columns([phenotype_col, p_col])
        .require_numeric([p_col])
        .require_range(p_col, min_val=0, max_val=1, exclusive_min=True)
        .validate()
    )

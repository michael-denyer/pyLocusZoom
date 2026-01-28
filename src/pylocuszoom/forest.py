"""Forest plot data validation and preparation.

Validates and prepares meta-analysis/forest plot data for visualization.
"""

import pandas as pd

from .validation import DataFrameValidator


def validate_forest_df(
    df: pd.DataFrame,
    study_col: str = "study",
    effect_col: str = "effect",
    ci_lower_col: str = "ci_lower",
    ci_upper_col: str = "ci_upper",
) -> None:
    """Validate forest plot DataFrame has required columns and types.

    Args:
        df: Forest plot data DataFrame.
        study_col: Column name for study/phenotype names.
        effect_col: Column name for effect sizes (beta, OR, HR).
        ci_lower_col: Column name for lower confidence interval.
        ci_upper_col: Column name for upper confidence interval.

    Raises:
        ValidationError: If required columns are missing or have invalid types.
    """
    (
        DataFrameValidator(df, "Forest plot DataFrame")
        .require_columns([study_col, effect_col, ci_lower_col, ci_upper_col])
        .require_numeric([effect_col, ci_lower_col, ci_upper_col])
        .require_ci_ordering(ci_lower_col, effect_col, ci_upper_col)
        .validate()
    )

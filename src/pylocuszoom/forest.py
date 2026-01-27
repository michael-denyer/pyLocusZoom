"""Forest plot data validation and preparation.

Validates and prepares meta-analysis/forest plot data for visualization.
"""

import pandas as pd

from .utils import ValidationError


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
    errors = []

    # Check required columns exist
    required = [study_col, effect_col, ci_lower_col, ci_upper_col]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValidationError(
            f"Forest plot DataFrame missing required columns: {missing}. "
            f"Required: {required}. Found: {list(df.columns)}"
        )

    # Check numeric columns are actually numeric
    numeric_cols = [effect_col, ci_lower_col, ci_upper_col]
    for col in numeric_cols:
        if not pd.api.types.is_numeric_dtype(df[col]):
            errors.append(f"Column '{col}' must be numeric, got {df[col].dtype}")

    if errors:
        raise ValidationError(
            "Forest plot validation failed:\n  - " + "\n  - ".join(errors)
        )

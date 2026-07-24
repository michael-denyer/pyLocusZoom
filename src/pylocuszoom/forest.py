"""Forest plot data validation and preparation.

Validates and prepares meta-analysis/forest plot data for visualization.
"""

import pandas as pd

from .exceptions import ForestValidationError
from .validation import ColumnSpec, check


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
        ForestValidationError: If required columns are missing or have invalid types.
    """
    check(
        df,
        ColumnSpec(
            name="Forest plot DataFrame",
            required=(study_col, effect_col, ci_lower_col, ci_upper_col),
            numeric=(effect_col, ci_lower_col, ci_upper_col),
            ordering=(
                (ci_lower_col, effect_col),
                (effect_col, ci_upper_col),
                (ci_lower_col, ci_upper_col),
            ),
            error_class=ForestValidationError,
        ),
    )

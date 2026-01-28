"""Tests for forest plot data validation."""

import pandas as pd
import pytest

from pylocuszoom.utils import ValidationError


def test_validate_forest_df_valid():
    """Test validation passes for valid forest plot DataFrame."""
    from pylocuszoom.forest import validate_forest_df

    df = pd.DataFrame(
        {
            "study": ["Study A", "Study B", "Study C"],
            "effect": [0.5, -0.2, 0.3],
            "ci_lower": [0.2, -0.5, 0.1],
            "ci_upper": [0.8, 0.1, 0.5],
        }
    )
    # Should not raise
    validate_forest_df(df)


def test_validate_forest_df_missing_column():
    """Test validation fails for missing required column."""
    from pylocuszoom.forest import validate_forest_df

    df = pd.DataFrame(
        {
            "study": ["Study A", "Study B"],
            "effect": [0.5, -0.2],
            # missing ci_lower and ci_upper
        }
    )
    with pytest.raises(ValidationError, match="ci_lower"):
        validate_forest_df(df)


def test_validate_forest_df_with_weight():
    """Test validation allows optional weight column."""
    from pylocuszoom.forest import validate_forest_df

    df = pd.DataFrame(
        {
            "study": ["Study A", "Study B"],
            "effect": [0.5, -0.2],
            "ci_lower": [0.2, -0.5],
            "ci_upper": [0.8, 0.1],
            "weight": [100, 200],
        }
    )
    # Should not raise
    validate_forest_df(df)


def test_validate_forest_df_ci_lower_gt_effect():
    """Test validation fails when ci_lower > effect.

    This would produce negative error bar lengths on the lower side.
    """
    from pylocuszoom.forest import validate_forest_df

    df = pd.DataFrame(
        {
            "study": ["Study A", "Study B"],
            "effect": [0.5, -0.2],
            "ci_lower": [0.7, -0.5],  # First row: 0.7 > 0.5 (effect)
            "ci_upper": [0.8, 0.1],
        }
    )
    with pytest.raises(ValidationError, match="ci_lower > effect"):
        validate_forest_df(df)


def test_validate_forest_df_effect_gt_ci_upper():
    """Test validation fails when effect > ci_upper.

    This would produce negative error bar lengths on the upper side.
    """
    from pylocuszoom.forest import validate_forest_df

    df = pd.DataFrame(
        {
            "study": ["Study A", "Study B"],
            "effect": [0.5, 0.3],
            "ci_lower": [0.2, 0.1],
            "ci_upper": [0.4, 0.5],  # First row: 0.5 > 0.4 (upper)
        }
    )
    with pytest.raises(ValidationError, match="effect > ci_upper"):
        validate_forest_df(df)


def test_validate_forest_df_ci_lower_gt_ci_upper():
    """Test validation fails when ci_lower > ci_upper.

    Completely inverted interval - no valid visualization.
    """
    from pylocuszoom.forest import validate_forest_df

    df = pd.DataFrame(
        {
            "study": ["Study A", "Study B"],
            "effect": [0.5, -0.2],
            "ci_lower": [0.8, -0.5],  # 0.8 > 0.2 (upper)
            "ci_upper": [0.2, 0.1],
        }
    )
    with pytest.raises(ValidationError, match="ci_lower > ci_upper"):
        validate_forest_df(df)

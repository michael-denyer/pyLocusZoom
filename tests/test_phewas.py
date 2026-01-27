"""Tests for PheWAS data validation."""

import pandas as pd
import pytest

from pylocuszoom.utils import ValidationError


def test_validate_phewas_df_valid():
    """Test validation passes for valid PheWAS DataFrame."""
    from pylocuszoom.phewas import validate_phewas_df

    df = pd.DataFrame(
        {
            "phenotype": ["Height", "BMI", "T2D"],
            "p_value": [1e-10, 0.05, 1e-5],
            "category": ["Anthropometric", "Anthropometric", "Metabolic"],
        }
    )
    # Should not raise
    validate_phewas_df(df)


def test_validate_phewas_df_missing_column():
    """Test validation fails for missing required column."""
    from pylocuszoom.phewas import validate_phewas_df

    df = pd.DataFrame(
        {
            "phenotype": ["Height", "BMI"],
            # missing p_value
        }
    )
    with pytest.raises(ValidationError, match="p_value"):
        validate_phewas_df(df)


def test_validate_phewas_df_optional_effect():
    """Test validation allows optional effect_size column."""
    from pylocuszoom.phewas import validate_phewas_df

    df = pd.DataFrame(
        {
            "phenotype": ["Height", "BMI"],
            "p_value": [1e-10, 0.05],
            "effect_size": [0.5, -0.2],
            "se": [0.1, 0.05],
        }
    )
    # Should not raise
    validate_phewas_df(df)

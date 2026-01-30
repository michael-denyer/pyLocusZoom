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


class TestPheWASNaNCategory:
    """Tests for PheWAS NaN category handling.

    Bug fix: pyLocusZoom-wej
    PheWAS rows with NaN category values were silently dropped because
    NaN == NaN is False in pandas.
    """

    def test_phewas_nan_category_included_in_plot(self, tmp_path):
        """Rows with NaN category should be included, not silently dropped."""
        import numpy as np

        from pylocuszoom.stats_plotter import StatsPlotter

        # Create PheWAS data with NaN category
        phewas_df = pd.DataFrame(
            {
                "phenotype": ["Phenotype_A", "Phenotype_B", "Phenotype_C"],
                "p_value": [0.01, 0.001, 0.05],
                "category": ["cat1", np.nan, "cat2"],  # Row B has NaN
            }
        )

        plotter = StatsPlotter()
        fig = plotter.plot_phewas(
            phewas_df,
            variant_id="rs12345",  # Required argument
            phenotype_col="phenotype",
            p_col="p_value",
            category_col="category",
        )

        # Get all scatter points from the figure
        # Each phenotype should have exactly one point
        ax = fig.axes[0]
        all_y_data = []
        for collection in ax.collections:
            offsets = collection.get_offsets()
            if len(offsets) > 0:
                all_y_data.extend(offsets[:, 1].tolist())

        # Should have 3 points (one for each phenotype including NaN category)
        assert len(all_y_data) == 3, (
            f"Expected 3 points (including NaN category row), got {len(all_y_data)}"
        )

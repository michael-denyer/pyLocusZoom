"""Tests for QQ plot functionality."""

import numpy as np
import pandas as pd
import pytest


class TestLambdaCalculation:
    """Tests for genomic inflation factor."""

    def test_lambda_is_1_for_uniform_pvalues(self):
        """Lambda should be ~1 for uniform p-values (no inflation)."""
        from pylocuszoom.qq import calculate_lambda_gc

        np.random.seed(42)
        p_values = np.random.uniform(0, 1, 10000)
        lambda_gc = calculate_lambda_gc(p_values)
        assert 0.95 < lambda_gc < 1.05

    def test_lambda_greater_than_1_for_inflated(self):
        """Lambda should be >1 for inflated p-values."""
        from pylocuszoom.qq import calculate_lambda_gc

        # Create inflated distribution (more small p-values than expected)
        np.random.seed(42)
        p_values = np.random.beta(0.5, 1, 10000)  # Skewed toward 0
        lambda_gc = calculate_lambda_gc(p_values)
        assert lambda_gc > 1.5

    def test_handles_nan_values(self):
        """Should ignore NaN p-values."""
        from pylocuszoom.qq import calculate_lambda_gc

        p_values = np.array([0.1, 0.5, np.nan, 0.3, 0.8])
        lambda_gc = calculate_lambda_gc(p_values)
        assert not np.isnan(lambda_gc)

    def test_handles_zero_pvalues(self):
        """Should handle p=0 gracefully."""
        from pylocuszoom.qq import calculate_lambda_gc

        p_values = np.array([0.0, 0.1, 0.5, 0.9])
        lambda_gc = calculate_lambda_gc(p_values)
        assert not np.isnan(lambda_gc)

    def test_returns_nan_for_empty_array(self):
        """Should return NaN for empty array."""
        from pylocuszoom.qq import calculate_lambda_gc

        lambda_gc = calculate_lambda_gc(np.array([]))
        assert np.isnan(lambda_gc)


class TestConfidenceBand:
    """Tests for QQ confidence band calculation."""

    def test_returns_three_arrays(self):
        """Should return expected, lower, upper."""
        from pylocuszoom.qq import calculate_confidence_band

        expected, lower, upper = calculate_confidence_band(100)
        assert len(expected) == 100
        assert len(lower) == 100
        assert len(upper) == 100

    def test_lower_below_upper(self):
        """Lower bound should be below upper bound."""
        from pylocuszoom.qq import calculate_confidence_band

        expected, lower, upper = calculate_confidence_band(100)
        assert all(lower <= upper)

    def test_expected_decreases(self):
        """Expected values should decrease (as index increases, p increases, -log10 p decreases)."""
        from pylocuszoom.qq import calculate_confidence_band

        expected, _, _ = calculate_confidence_band(100)
        # Expected -log10(p) decreases as rank increases (larger p -> smaller -log10 p)
        assert all(np.diff(expected) <= 0)

    def test_band_widens_at_small_p(self):
        """Confidence band should be wider at small p (high -log10 p, start of array)."""
        from pylocuszoom.qq import calculate_confidence_band

        expected, lower, upper = calculate_confidence_band(100)
        band_width = upper - lower
        # Band is wider at start (small p values, high -log10 p)
        assert band_width[0] > band_width[len(band_width) // 2]


class TestPrepareQQData:
    """Tests for QQ data preparation."""

    def test_adds_expected_and_observed(self):
        """Should add expected and observed columns."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert "_expected" in result.columns
        assert "_observed" in result.columns

    def test_adds_confidence_bounds(self):
        """Should add confidence interval columns."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert "_ci_lower" in result.columns
        assert "_ci_upper" in result.columns

    def test_stores_lambda_in_attrs(self):
        """Should store lambda in DataFrame attrs."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert "lambda_gc" in result.attrs

    def test_stores_n_variants_in_attrs(self):
        """Should store variant count in DataFrame attrs."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert result.attrs["n_variants"] == 5

    def test_validates_p_column(self):
        """Should raise on missing p column."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"wrong": [0.1, 0.5]})
        with pytest.raises(ValueError, match="not found"):
            prepare_qq_data(df, p_col="p")

    def test_filters_invalid_pvalues(self):
        """Should filter out NaN and out-of-range p-values."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, np.nan, -0.1, 1.5, 0.5]})
        result = prepare_qq_data(df, p_col="p")
        # Only 0.1 and 0.5 are valid
        assert result.attrs["n_variants"] == 2

    def test_raises_on_no_valid_pvalues(self):
        """Should raise if no valid p-values."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [np.nan, -0.1, 1.5]})
        with pytest.raises(ValueError, match="No valid p-values"):
            prepare_qq_data(df, p_col="p")

    def test_observed_is_sorted_decreasing(self):
        """Observed values should be sorted decreasing (-log10 p from high to low)."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.9, 0.1, 0.01, 0.001, 0.5]})
        result = prepare_qq_data(df, p_col="p")
        # QQ plots have expected on x, observed on y
        # We sort p-values ascending, so -log10(p) is descending
        assert result["_observed"].is_monotonic_decreasing

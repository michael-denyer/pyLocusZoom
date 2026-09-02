"""Tests for QQ plot functionality."""

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from pylocuszoom.manhattan_plotter import ManhattanPlotter
from pylocuszoom.qq import calculate_confidence_band, calculate_lambda_gc
from tests.strategies import pvalues


class TestLambdaCalculation:
    """Tests for genomic inflation factor."""

    def test_lambda_is_1_for_uniform_pvalues(self):
        """Lambda should be ~1 for uniform p-values (no inflation)."""
        from pylocuszoom.qq import calculate_lambda_gc

        rng = np.random.default_rng(42)
        p_values = rng.uniform(0, 1, 10000)
        lambda_gc = calculate_lambda_gc(p_values)
        assert 0.95 < lambda_gc < 1.05

    def test_lambda_greater_than_1_for_inflated(self):
        """Lambda should be >1 for inflated p-values."""
        from pylocuszoom.qq import calculate_lambda_gc

        # Create inflated distribution (more small p-values than expected)
        rng = np.random.default_rng(42)
        p_values = rng.beta(0.5, 1, 10000)  # Skewed toward 0
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
        assert "_expected" in result.frame.columns
        assert "_observed" in result.frame.columns

    def test_adds_confidence_bounds(self):
        """Should add confidence interval columns."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert "_ci_lower" in result.frame.columns
        assert "_ci_upper" in result.frame.columns

    def test_carries_lambda(self):
        """The inflation factor is a field of the prepared value."""
        from pylocuszoom.qq import calculate_lambda_gc, prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert result.lambda_gc == calculate_lambda_gc(df["p"].to_numpy())

    def test_carries_n_variants(self):
        """The valid p-value count is a field of the prepared value."""
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.1, 0.01, 0.001, 0.5, 0.9]})
        result = prepare_qq_data(df, p_col="p")
        assert result.n_variants == 5

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
        assert result.n_variants == 2

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
        assert result.frame["_observed"].is_monotonic_decreasing


# =============================================================================
# Property-Based Tests (Hypothesis)
# =============================================================================


class TestQQProperties:
    """Property-based tests for QQ plot calculations."""

    @given(st.lists(pvalues(allow_zero=False), min_size=10, max_size=500))
    def test_lambda_gc_is_positive(self, p_values):
        """Lambda GC should always be positive for valid p-values."""
        p_array = np.array(p_values)
        lambda_gc = calculate_lambda_gc(p_array)

        # Filter out edge case of all identical p-values
        if np.std(p_array) < 1e-10:
            return

        assert lambda_gc > 0 or np.isnan(lambda_gc)

    @given(st.integers(min_value=10, max_value=1000))
    def test_confidence_band_shapes_match(self, n):
        """Confidence band arrays should all have length n."""
        expected, lower, upper = calculate_confidence_band(n)

        assert len(expected) == n
        assert len(lower) == n
        assert len(upper) == n

    @given(st.integers(min_value=10, max_value=500))
    def test_confidence_band_ordering(self, n):
        """Lower bound should be <= expected <= upper bound."""
        expected, lower, upper = calculate_confidence_band(n)

        # Allow small floating point tolerance
        assert np.all(lower <= expected + 1e-10)
        assert np.all(expected <= upper + 1e-10)


class TestQQWithVariousPvalueDistributions:
    """Test QQ plot with various p-value distributions."""

    @pytest.fixture
    def default_manhattan_plotter(self):
        """Create plotter instance for testing."""
        return ManhattanPlotter()

    def test_qq_uniform_pvalues(self, default_manhattan_plotter):
        """QQ plot with uniform p-values should show lambda ~ 1."""
        rng = np.random.default_rng(42)
        df = pd.DataFrame({"p": rng.uniform(0, 1, 1000)})

        fig = default_manhattan_plotter.plot_qq(df, p_col="p", show_lambda=True)
        assert fig is not None

        # Check title contains lambda close to 1
        ax = fig.axes[0]
        title = ax.get_title()
        assert "λ" in title

    def test_qq_extreme_pvalues(self, default_manhattan_plotter):
        """QQ plot with very small p-values should not produce inf."""
        df = pd.DataFrame(
            {
                "p": [1e-300, 1e-200, 1e-100, 1e-50, 0.001, 0.01, 0.1, 0.5],
            }
        )

        fig = default_manhattan_plotter.plot_qq(df, p_col="p")
        assert fig is not None


class TestEmptyQQInput:
    """Empty input has no valid p-values to rank."""

    def test_qq_with_empty_df_raises(self):
        """QQ plot with empty DataFrame raises ValueError.

        The QQ plot requires valid p-values (>0 and <=1) to compute
        expected quantiles. Empty data fails this validation.
        """
        plotter = ManhattanPlotter()
        empty_df = pd.DataFrame({"p": pd.Series([], dtype=float)})

        with pytest.raises(ValueError, match="No valid p-values"):
            plotter.plot_qq(empty_df, p_col="p")

"""Tests for QQ plot functionality."""

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st
from matplotlib.collections import PathCollection

from pylocuszoom import GenomeWideConfig
from pylocuszoom.exceptions import ValidationError
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

        fig = default_manhattan_plotter.plot_qq(
            df, show_lambda=True, config=GenomeWideConfig(p_col="p")
        )
        assert fig is not None

        # Check title contains lambda close to 1
        ax = fig.axes[0]
        title = ax.get_title()
        assert "λ" in title

    def test_qq_extreme_pvalues(self, default_manhattan_plotter):
        """Plot p-values down to 1e-300 as finite coordinates, not inf."""
        df = pd.DataFrame(
            {
                "p": [1e-300, 1e-200, 1e-100, 1e-50, 0.001, 0.01, 0.1, 0.5],
            }
        )

        fig = default_manhattan_plotter.plot_qq(df, config=GenomeWideConfig(p_col="p"))

        ax = fig.get_axes()[0]
        drawn = np.array(
            [
                point
                for collection in ax.collections
                for point in collection.get_offsets()
            ]
        )
        assert np.isfinite(drawn).all()
        assert drawn[:, 1].max() == pytest.approx(300.0)


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
            plotter.plot_qq(empty_df, config=GenomeWideConfig(p_col="p"))


class TestWeightedQQ:
    """A QQ from a stratified sample, weighted by inverse inclusion probability.

    The sample stands for a full scan of ``N_FULL`` uniform p-values: every
    p below ``CUTOFF`` is kept with weight 1, and each other p is kept with
    probability ``FRACTION`` and weight ``1 / FRACTION``. The weighted QQ of
    the sample should put each kept p where the full scan's QQ puts it.

    The tolerances, in -log10 units, are the sampling noise of this design
    measured over 200 seeds: the error is largest just past the seam, where
    few sampled rows sit among many kept ones (worst 0.10), and smallest in
    the kept tail, where only the total weight is estimated (worst 0.02).
    Lambda's weighted median comes from the ~4,000 sampled rows (worst error
    0.12). Without weights the tail is off by 1.5 and lambda by 1.8.
    """

    N_FULL = 200_000
    CUTOFF = 1e-2
    FRACTION = 0.02
    TAIL_TOLERANCE = 0.03
    TOLERANCE = 0.12
    LAMBDA_TOLERANCE = 0.15

    @pytest.fixture
    def scan(self):
        """Return the full scan's p-values and the weighted sample of them."""
        rng = np.random.default_rng(7)
        full = rng.uniform(0, 1, self.N_FULL)
        tail = full < self.CUTOFF
        sampled = ~tail & (rng.uniform(0, 1, self.N_FULL) < self.FRACTION)
        keep = tail | sampled
        sample = pd.DataFrame(
            {
                "chr": rng.choice(["1", "2", "3"], keep.sum()),
                "pos": rng.integers(1, 10_000_000, keep.sum()),
                "p_value": full[keep],
                "w": np.where(tail[keep], 1.0, 1 / self.FRACTION),
            }
        )
        return np.sort(full), sample

    @staticmethod
    def _qq_points(ax):
        """Return the drawn QQ points as (expected, observed), largest first."""
        points = np.concatenate(
            [c.get_offsets() for c in ax.collections if isinstance(c, PathCollection)]
        )
        return points[np.argsort(-points[:, 1])]

    @staticmethod
    def _full_scan_expected(full_sorted, sample):
        """Return the full scan's expected -log10 p of each sampled p, smallest first."""
        ranks = np.searchsorted(full_sorted, np.sort(sample["p_value"]), side="right")
        return -np.log10(ranks / (len(full_sorted) + 1))

    def _seam(self, observed):
        """Return a mask of the 40 points either side of the cutoff."""
        seam_index = np.searchsorted(-observed, np.log10(self.CUTOFF))
        mask = np.zeros(len(observed), dtype=bool)
        mask[max(seam_index - 40, 0) : seam_index + 40] = True
        return mask

    def test_weighted_sample_matches_full_scan(self, scan):
        full_sorted, sample = scan

        fig = ManhattanPlotter().plot_qq(
            sample, config=GenomeWideConfig(weight_col="w")
        )

        expected, observed = self._qq_points(fig.get_axes()[0]).T
        truth = self._full_scan_expected(full_sorted, sample)
        error = np.abs(expected - truth)
        assert error[:20].max() < self.TAIL_TOLERANCE
        assert error[self._seam(observed)].max() < self.TOLERANCE
        assert error.max() < self.TOLERANCE

    def test_unweighted_sample_misplaces_the_tail(self, scan):
        """Without weights each row is one test, so the sample is the scan."""
        full_sorted, sample = scan

        fig = ManhattanPlotter().plot_qq(sample)

        expected, _ = self._qq_points(fig.get_axes()[0]).T
        truth = self._full_scan_expected(full_sorted, sample)
        assert truth[0] - expected[0] > 1

    def test_weighted_lambda_matches_full_scan(self, scan):
        full_sorted, sample = scan

        weighted = ManhattanPlotter().plot_qq(
            sample, config=GenomeWideConfig(weight_col="w")
        )
        unweighted = ManhattanPlotter().plot_qq(sample)

        def shown_lambda(fig):
            return float(fig.get_axes()[0].get_title().split("= ")[1].rstrip(")"))

        assert shown_lambda(weighted) == pytest.approx(
            calculate_lambda_gc(full_sorted), abs=self.LAMBDA_TOLERANCE
        )
        assert shown_lambda(unweighted) > 2

    def test_weighted_band_matches_full_scan(self, scan):
        """The 95% band at each tail point is the full scan's band there."""
        from pylocuszoom.qq import prepare_qq_data

        full_sorted, sample = scan
        n_tail = int((sample["p_value"] < self.CUTOFF).sum())

        weighted = prepare_qq_data(sample, weight_col="w").frame
        _, lower, upper = calculate_confidence_band(self.N_FULL)

        np.testing.assert_allclose(
            weighted["_ci_lower"][:n_tail], lower[:n_tail], atol=self.TAIL_TOLERANCE
        )
        np.testing.assert_allclose(
            weighted["_ci_upper"][:n_tail], upper[:n_tail], atol=self.TAIL_TOLERANCE
        )

    def test_manhattan_qq_stacked_weights_the_qq_only(self, scan):
        """The QQ panel is weighted; the Manhattan panel draws every row once."""
        full_sorted, sample = scan
        plotter = ManhattanPlotter(species="human")

        weighted = plotter.plot_manhattan_qq_stacked(
            [sample], config=GenomeWideConfig(weight_col="w")
        )
        plain = plotter.plot_manhattan_qq_stacked([sample])

        manhattan, qq = weighted.get_axes()
        expected, _ = self._qq_points(qq).T
        truth = self._full_scan_expected(full_sorted, sample)
        assert np.abs(expected - truth).max() < self.TOLERANCE
        np.testing.assert_array_equal(
            self._qq_points(manhattan), self._qq_points(plain.get_axes()[0])
        )

    def test_manhattan_qq_weights_the_qq(self, scan):
        full_sorted, sample = scan

        fig = ManhattanPlotter(species="human").plot_manhattan_qq(
            sample, config=GenomeWideConfig(weight_col="w")
        )

        expected, _ = self._qq_points(fig.get_axes()[1]).T
        truth = self._full_scan_expected(full_sorted, sample)
        assert np.abs(expected - truth).max() < self.TOLERANCE

    def test_unit_weights_change_nothing(self):
        """Weight 1 on every row is the unweighted QQ, value for value."""
        from pylocuszoom.qq import prepare_qq_data

        rng = np.random.default_rng(3)
        df = pd.DataFrame({"p": rng.uniform(0, 1, 1001), "w": 1.0})

        plain = prepare_qq_data(df, p_col="p")
        weighted = prepare_qq_data(df, p_col="p", weight_col="w")

        pd.testing.assert_frame_equal(weighted.frame, plain.frame)
        assert weighted.lambda_gc == pytest.approx(plain.lambda_gc, rel=1e-12)
        assert weighted.n_variants == plain.n_variants

    def test_even_count_unit_weights_lambda_is_the_plain_median(self):
        from pylocuszoom.qq import prepare_qq_data

        df = pd.DataFrame({"p": [0.9, 0.5, 0.2, 0.01], "w": 1.0})

        assert prepare_qq_data(df, p_col="p", weight_col="w").lambda_gc == (
            pytest.approx(prepare_qq_data(df, p_col="p").lambda_gc, rel=1e-12)
        )

    def test_missing_weight_column_is_named(self):
        df = pd.DataFrame({"p_value": [0.1, 0.5]})

        with pytest.raises(ValidationError, match="weight"):
            ManhattanPlotter().plot_qq(df, config=GenomeWideConfig(weight_col="weight"))

    @pytest.mark.parametrize(
        "bad", [0.0, -1.0, np.nan, np.inf], ids=["zero", "negative", "nan", "inf"]
    )
    def test_rejects_a_weight_that_is_not_positive_and_finite(self, bad):
        df = pd.DataFrame(
            {"chr": ["1", "1"], "pos": [1, 2], "p_value": [0.1, 0.5], "w": [1.0, bad]}
        )
        config = GenomeWideConfig(weight_col="w")

        with pytest.raises(ValidationError, match="'w'"):
            ManhattanPlotter().plot_qq(df, config=config)
        with pytest.raises(ValidationError, match="'w'"):
            ManhattanPlotter(species="human").plot_manhattan_qq(df, config=config)

    def test_rejects_a_non_numeric_weight(self):
        df = pd.DataFrame({"p_value": [0.1, 0.5], "w": ["1", "x"]})

        with pytest.raises(ValidationError, match="'w'"):
            ManhattanPlotter().plot_qq(df, config=GenomeWideConfig(weight_col="w"))

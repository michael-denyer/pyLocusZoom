"""QQ plot data preparation and statistics."""

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

from ._data import prepare_pvalue_data
from .exceptions import ValidationError
from .schemas import Canonical, qq_weight_spec
from .validation import check


def calculate_lambda_gc(
    p_values: np.ndarray, weights: Optional[np.ndarray] = None
) -> float:
    """Calculate genomic inflation factor (lambda GC).

    Lambda is the ratio of the median observed chi-squared statistic
    to the expected median under the null hypothesis.

    Args:
        p_values: Array of p-values.
        weights: Number of tests each p-value stands for, aligned with
            ``p_values``, or None for one each. The median is then the
            weighted median, so a stratified sample gives the lambda of the
            scan it was drawn from.

    Returns:
        Genomic inflation factor (lambda). Returns NaN if no valid p-values.
    """
    # Remove NaN and zero/negative values
    valid = ~np.isnan(p_values) & (p_values > 0)
    p_clean = p_values[valid]
    if len(p_clean) == 0:
        return np.nan

    # Convert to chi-squared statistics (1 df)
    chi2 = stats.chi2.ppf(1 - p_clean, df=1)

    # Expected median for chi-squared with 1 df
    expected_median = stats.chi2.ppf(0.5, df=1)

    if weights is None:
        observed_median = np.median(chi2)
    else:
        observed_median = _weighted_median(chi2, weights[valid])

    # Lambda = observed median / expected median
    return observed_median / expected_median


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Return the weighted median, which is ``np.median`` for equal weights."""
    order = np.argsort(values)
    values, cumulative = values[order], np.cumsum(weights[order])
    half = cumulative[-1] / 2
    index = int(np.searchsorted(cumulative, half))
    if np.isclose(cumulative[index], half) and index + 1 < len(values):
        return (values[index] + values[index + 1]) / 2
    return values[index]


def calculate_confidence_band(
    n_points: int, confidence: float = 0.95
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculate confidence band for QQ plot.

    Uses order statistics to compute expected distribution of p-values
    under the null hypothesis.

    Args:
        n_points: Number of p-values.
        confidence: Confidence level (default 0.95 for 95% CI).

    Returns:
        Tuple of (expected, lower_bound, upper_bound) arrays in -log10 scale.
    """
    return _rank_quantiles(np.arange(1, n_points + 1), n_points, confidence)


def _rank_quantiles(
    ranks: np.ndarray, n_tests: float, confidence: float = 0.95
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the expected -log10 p and its band at each rank among ``n_tests``.

    Ranks and ``n_tests`` may be fractional: ``scipy.stats.beta`` takes real
    shape parameters, and a weighted rank is an estimate of a whole one.
    """
    # Expected quantiles
    expected = -np.log10(ranks / (n_tests + 1))

    # Confidence interval using beta distribution
    alpha = 1 - confidence
    n_minus_rank = n_tests - ranks + 1

    lower_p = stats.beta.ppf(alpha / 2, ranks, n_minus_rank)
    upper_p = stats.beta.ppf(1 - alpha / 2, ranks, n_minus_rank)

    # Convert to -log10 scale (swap because -log10 reverses order)
    lower_bound = -np.log10(upper_p)
    upper_bound = -np.log10(lower_p)

    return expected, lower_bound, upper_bound


@dataclass(frozen=True)
class PreparedQQ:
    """Observed and expected quantiles for one QQ panel, with its statistics.

    ``frame`` carries ``_expected``, ``_observed``, ``_ci_lower`` and
    ``_ci_upper``; ``lambda_gc`` is the genomic inflation factor and
    ``n_variants`` the number of valid p-values behind them.
    """

    frame: pd.DataFrame
    lambda_gc: float
    n_variants: int


def prepare_qq_data(
    df: pd.DataFrame,
    p_col: str = Canonical.P,
    weight_col: Optional[str] = None,
) -> PreparedQQ:
    """Prepare DataFrame for QQ plot rendering.

    Without weights each row is one test: the ``i``-th smallest of ``n``
    p-values is expected at ``-log10(i / (n + 1))``. With weights each row
    stands for ``weight`` tests, so a stratified sample (every p below a
    cutoff with weight 1, a random fraction ``f`` of the rest with weight
    ``1 / f``) is drawn where the whole scan would be. A row's rank is the
    cumulative weight up to and including it in p order, and ``n`` is the
    total weight, so unit weights give ``i / (n + 1)`` exactly. A kept row's
    rank is then its rank in the scan, and a sampled row's is an unbiased
    estimate of it: the ``k``-th row kept at random with probability ``f``
    sits ``k / f`` ranks into the rows it was drawn from, on average. The
    band and lambda use the same ranks and weights.

    Args:
        df: DataFrame with p-values.
        p_col: Column name for p-value.
        weight_col: Column holding the number of tests each row stands for,
            or None for one each.

    Returns:
        The quantile frame and its statistics. ``n_variants`` counts rows,
        not weight.

    Raises:
        ValidationError: If ``p_col`` is missing, if no p-value lies in
            ``(0, 1]``, or if ``weight_col`` is missing or holds a value
            that is null, non-numeric, or not positive and finite.
    """
    if p_col not in df.columns:
        raise ValidationError(f"Column '{p_col}' not found in DataFrame")
    if weight_col is not None:
        check(df, qq_weight_spec(weight_col))

    prepared = prepare_pvalue_data(
        df,
        p_col,
        "qq",
        on_empty="No valid p-values found (must be > 0 and <= 1)",
    )
    p_valid = pd.to_numeric(prepared[p_col]).to_numpy()

    # Sort p-values (smallest first -> largest -log10 last)
    order = np.argsort(p_valid, kind="stable")
    p_sorted = p_valid[order]

    # Calculate observed -log10(p)
    observed = -np.log10(p_sorted)

    # Calculate expected and confidence bands
    if weight_col is None:
        weights = None
        expected, ci_lower, ci_upper = calculate_confidence_band(len(p_sorted))
    else:
        weights = prepared[weight_col].to_numpy(dtype=float)
        sorted_weights = weights[order]
        expected, ci_lower, ci_upper = _rank_quantiles(
            np.cumsum(sorted_weights), sorted_weights.sum()
        )

    frame = pd.DataFrame(
        {
            "_expected": expected,
            "_observed": observed,
            "_ci_lower": ci_lower,
            "_ci_upper": ci_upper,
        }
    )
    return PreparedQQ(frame, calculate_lambda_gc(p_valid, weights), len(p_valid))

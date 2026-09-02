"""QQ plot data preparation and statistics."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from ._data import prepare_pvalue_data


def calculate_lambda_gc(p_values: np.ndarray) -> float:
    """Calculate genomic inflation factor (lambda GC).

    Lambda is the ratio of the median observed chi-squared statistic
    to the expected median under the null hypothesis.

    Args:
        p_values: Array of p-values.

    Returns:
        Genomic inflation factor (lambda). Returns NaN if no valid p-values.
    """
    # Remove NaN and zero/negative values
    p_clean = p_values[~np.isnan(p_values) & (p_values > 0)]
    if len(p_clean) == 0:
        return np.nan

    # Convert to chi-squared statistics (1 df)
    chi2 = stats.chi2.ppf(1 - p_clean, df=1)

    # Expected median for chi-squared with 1 df
    expected_median = stats.chi2.ppf(0.5, df=1)

    # Lambda = observed median / expected median
    return np.median(chi2) / expected_median


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
    # Expected quantiles
    expected = -np.log10((np.arange(1, n_points + 1)) / (n_points + 1))

    # Confidence interval using beta distribution
    alpha = 1 - confidence
    ranks = np.arange(1, n_points + 1)
    n_minus_rank = n_points - ranks + 1

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
    p_col: str = "p",
) -> PreparedQQ:
    """Prepare DataFrame for QQ plot rendering.

    Args:
        df: DataFrame with p-values.
        p_col: Column name for p-value.

    Returns:
        The quantile frame and its statistics.

    Raises:
        ValueError: If ``p_col`` is missing or no p-value lies in ``(0, 1]``.
    """
    if p_col not in df.columns:
        raise ValueError(f"Column '{p_col}' not found in DataFrame")

    prepared = prepare_pvalue_data(
        df,
        p_col,
        allow_zero=False,
        on_empty="No valid p-values found (must be > 0 and <= 1)",
    )
    p_valid = pd.to_numeric(prepared[p_col]).to_numpy()

    # Sort p-values (smallest first -> largest -log10 last)
    p_sorted = np.sort(p_valid)

    # Calculate observed -log10(p)
    observed = -np.log10(p_sorted)

    # Calculate expected and confidence bands
    expected, ci_lower, ci_upper = calculate_confidence_band(len(p_sorted))

    frame = pd.DataFrame(
        {
            "_expected": expected,
            "_observed": observed,
            "_ci_lower": ci_lower,
            "_ci_upper": ci_upper,
        }
    )
    return PreparedQQ(frame, calculate_lambda_gc(p_valid), len(p_valid))

"""QQ plot data preparation and statistics."""

import numpy as np
import pandas as pd
from scipy import stats


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


def prepare_qq_data(
    df: pd.DataFrame,
    p_col: str = "p",
) -> pd.DataFrame:
    """Prepare DataFrame for QQ plot rendering.

    Args:
        df: DataFrame with p-values.
        p_col: Column name for p-value.

    Returns:
        DataFrame with columns for QQ plotting:
        - _expected: Expected -log10(p) under null
        - _observed: Observed -log10(p)
        - _ci_lower: Lower confidence bound
        - _ci_upper: Upper confidence bound

        Attributes stored in DataFrame.attrs:
        - lambda_gc: Genomic inflation factor
        - n_variants: Number of valid p-values
    """
    if p_col not in df.columns:
        raise ValueError(f"Column '{p_col}' not found in DataFrame")

    # Get p-values and filter invalid
    p_values = df[p_col].values
    valid_mask = ~np.isnan(p_values) & (p_values > 0) & (p_values <= 1)
    p_valid = p_values[valid_mask]

    if len(p_valid) == 0:
        raise ValueError("No valid p-values found (must be > 0 and <= 1)")

    # Sort p-values (smallest first -> largest -log10 last)
    p_sorted = np.sort(p_valid)

    # Calculate observed -log10(p)
    observed = -np.log10(p_sorted)

    # Calculate expected and confidence bands
    expected, ci_lower, ci_upper = calculate_confidence_band(len(p_sorted))

    # Create result DataFrame
    result = pd.DataFrame(
        {
            "_expected": expected,
            "_observed": observed,
            "_ci_lower": ci_lower,
            "_ci_upper": ci_upper,
        }
    )

    # Store statistics in attrs
    result.attrs["lambda_gc"] = calculate_lambda_gc(p_valid)
    result.attrs["n_variants"] = len(p_valid)

    return result

"""Shared utilities for plotter classes.

Internal module - not part of public API.
"""

from typing import Any, Optional

import numpy as np
import pandas as pd

# Significance thresholds
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8

# Manhattan/QQ plot styling constants
MANHATTAN_POINT_SIZE = 10
MANHATTAN_CATEGORICAL_POINT_SIZE = 30
QQ_POINT_SIZE = 10
POINT_EDGE_COLOR = "black"
MANHATTAN_EDGE_WIDTH = 0.1
QQ_EDGE_WIDTH = 0.02
QQ_POINT_COLOR = "#1f77b4"
QQ_CI_COLOR = "#CCCCCC"
QQ_CI_ALPHA = 0.5
SIGNIFICANCE_LINE_COLOR = "red"


def transform_pvalues(df: pd.DataFrame, p_col: str) -> pd.DataFrame:
    """Validate, filter, and -log10 transform p-values.

    Filters out invalid rows before transformation:
    - NaN p-values are removed.
    - Out-of-range p-values (< 0 or > 1) are removed.
    - Very small valid p-values (< 1e-300) are clipped to avoid -inf.

    Args:
        df: DataFrame with p-value column.
        p_col: Name of p-value column.

    Returns:
        DataFrame with invalid rows removed and neglog10p column added.
    """
    from .logging import logger

    df = df.copy()

    # Filter NaN p-values
    nan_mask = df[p_col].isna()
    nan_count = nan_mask.sum()
    if nan_count > 0:
        logger.warning("Found {} NaN p-values, filtering out", nan_count)
        df = df[~nan_mask]

    # Filter out-of-range p-values
    invalid_mask = (df[p_col] < 0) | (df[p_col] > 1)
    invalid_count = invalid_mask.sum()
    if invalid_count > 0:
        logger.warning(
            "Found {} p-values outside [0, 1] range, filtering out",
            invalid_count,
        )
        df = df[~invalid_mask]

    df["neglog10p"] = -np.log10(df[p_col].clip(lower=1e-300))
    return df


def add_significance_line(
    backend: Any,
    ax: Any,
    threshold: Optional[float],
) -> None:
    """Add genome-wide significance threshold line.

    Args:
        backend: Plot backend instance.
        ax: Axes object from backend.
        threshold: P-value threshold (e.g., 5e-8). None to skip.
    """
    if threshold is None:
        return
    threshold_line = -np.log10(threshold)
    backend.axhline(
        ax,
        y=threshold_line,
        color=SIGNIFICANCE_LINE_COLOR,
        linestyle="--",
        linewidth=1,
        zorder=1,
    )

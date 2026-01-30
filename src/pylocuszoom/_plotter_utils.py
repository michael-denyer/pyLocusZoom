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
    """Add neglog10p column with -log10 transformed p-values.

    Clips extremely small p-values to 1e-300 to avoid -inf.

    Args:
        df: DataFrame with p-value column.
        p_col: Name of p-value column.

    Returns:
        DataFrame with neglog10p column added.
    """
    df = df.copy()
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

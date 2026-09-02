"""Shared plot-data intake policy."""

from typing import Optional

import numpy as np
import pandas as pd

from .logging import logger

P_VALUE_FLOOR = 1e-300
# Upper bound of the p-value domain, shared with ColumnSpec.pvalue so strict
# validation and plot-time intake cannot drift apart.
P_VALUE_MAX = 1.0


def prepare_pvalue_data(
    df: pd.DataFrame,
    p_col: str,
    *,
    allow_zero: bool = True,
    out_col: str = "neglog10p",
    on_empty: Optional[str] = None,
) -> pd.DataFrame:
    """Return a copy with valid p-values and a finite ``-log10(p)`` column.

    Null, non-numeric, out-of-range p-values are filtered consistently for
    every plot family. Tiny valid values are clipped before taking ``-log10``.

    Args:
        df: Frame holding the p-value column.
        p_col: Name of the p-value column.
        allow_zero: Keep an exact zero as a valid, clipped p-value, the
            Manhattan convention. Strict callers such as eQTL and QQ opt into
            the mathematical ``(0, 1]`` domain by passing False.
        out_col: Name of the transformed column to write.
        on_empty: Message to raise when no row survives filtering. None
            returns the empty frame instead, which the regional path relies on.

    Returns:
        A filtered copy carrying ``out_col``.

    Raises:
        ValueError: If nothing survives and ``on_empty`` names a message.
    """
    result = df.copy()
    initial_count = len(result)
    p_values = pd.to_numeric(result[p_col], errors="coerce")
    lower_mask = p_values >= 0 if allow_zero else p_values > 0
    valid = p_values.notna() & lower_mask & (p_values <= P_VALUE_MAX)
    dropped = int((~valid).sum())
    nan_count = int(result[p_col].isna().sum())
    if nan_count:
        logger.warning("Found {} NaN p-values, filtering out", nan_count)
    numeric_missing = int(p_values.isna().sum()) - nan_count
    if numeric_missing:
        logger.warning("Found {} non-numeric p-values, filtering out", numeric_missing)
    out_of_range = int((p_values.notna() & ~valid).sum())
    if out_of_range:
        logger.warning(
            "Found {} p-values outside [0, 1] range, filtering out", out_of_range
        )
    result = result.loc[valid].copy()
    valid_values = p_values.loc[valid]
    clipped = int((valid_values < P_VALUE_FLOOR).sum())
    if clipped:
        logger.debug("Clipping {} p-values below {}", clipped, P_VALUE_FLOOR)
    result[out_col] = -np.log10(valid_values.clip(lower=P_VALUE_FLOOR))
    if dropped:
        logger.debug("P-value filtering removed {} of {} rows", dropped, initial_count)
    if result.empty and on_empty is not None:
        raise ValueError(on_empty)
    return result

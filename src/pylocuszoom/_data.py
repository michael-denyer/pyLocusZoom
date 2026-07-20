"""Shared plot-data intake policy."""

import numpy as np
import pandas as pd

from .logging import logger

P_VALUE_FLOOR = 1e-300


def prepare_pvalue_data(
    df: pd.DataFrame, p_col: str, *, allow_zero: bool = True
) -> pd.DataFrame:
    """Return a copy with valid p-values and a finite ``neglog10p`` column.

    Null, non-numeric, out-of-range p-values are filtered consistently for
    every plot family. ``allow_zero`` preserves the Manhattan convention that
    an exact zero is a valid, clipped p-value; strict callers such as eQTL can
    opt into the mathematical ``(0, 1]`` domain. Tiny valid values are clipped
    before taking ``-log10``.
    """
    result = df.copy()
    initial_count = len(result)
    p_values = pd.to_numeric(result[p_col], errors="coerce")
    lower_bound = 0 if allow_zero else 0
    lower_mask = p_values >= lower_bound if allow_zero else p_values > lower_bound
    valid = p_values.notna() & lower_mask & (p_values <= 1)
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
    result["neglog10p"] = -np.log10(valid_values.clip(lower=P_VALUE_FLOOR))
    if dropped:
        logger.debug("P-value filtering removed {} of {} rows", dropped, initial_count)
    return result

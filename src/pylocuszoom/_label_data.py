"""Backend-neutral selection of SNP rows eligible for text labels."""

import pandas as pd

from .exceptions import ValidationError


def select_label_candidates(
    df: pd.DataFrame,
    *,
    pos_col: str,
    region_span: int | None,
    lead_pos: int | None = None,
    lead_index: int | None = None,
    min_label_distance: float = 0.05,
) -> pd.DataFrame:
    """Keep the lead and rows outside its proximity exclusion zone.

    Regional callers identify one row with ``lead_index``. The standalone
    label API supplies ``lead_pos`` instead, preserving its contract that
    every row at that position is exempt. An index takes precedence and
    supplies the position from that row. Without a lead or positive region
    span, every row remains eligible. Selection precedes top-N ranking.
    """
    if lead_index is not None:
        lead_pos = df.at[lead_index, pos_col]
    if lead_pos is None or region_span is None or region_span <= 0:
        return df
    if not 0 <= min_label_distance <= 1:
        raise ValidationError(
            f"min_label_distance must be between 0 and 1, got {min_label_distance}"
        )
    is_lead = df[pos_col].eq(lead_pos) if lead_index is None else df.index == lead_index
    distance = (df[pos_col] - lead_pos).abs()
    return df[is_lead | distance.ge(min_label_distance * region_span)]

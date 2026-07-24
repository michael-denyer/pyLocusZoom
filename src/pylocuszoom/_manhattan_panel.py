"""Shared Manhattan-panel drawing primitives.

A Miami plot is a mirrored Manhattan, so the single/stacked renderer
(:mod:`._rendering`) and the Miami renderer (:mod:`._miami_renderer`) share the
same per-chromosome scatter loop, y-limit padding, and cumulative-position
x-limit padding. Those primitives live here so both renderers stay in step.
"""

from typing import Any, Optional, Sequence, Tuple

import pandas as pd

from ._plotter_utils import (
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
)
from .backends.base import PlotBackend
from .backends.hover import HoverConfig, HoverDataBuilder


def padded_ymax(y_max: float) -> float:
    """Return a useful upper y-limit for a Manhattan panel."""
    return max(y_max * 1.1, 1.0) if pd.notna(y_max) else 1.0


def shared_manhattan_limits(
    prepared_dfs: Sequence[pd.DataFrame],
) -> Tuple[float, float]:
    """Return padded x-limits spanning the cumulative positions of every panel."""
    x_min = min(df["_cumulative_pos"].min() for df in prepared_dfs)
    x_max = max(df["_cumulative_pos"].max() for df in prepared_dfs)
    x_padding = (x_max - x_min) * 0.01
    return x_min - x_padding, x_max + x_padding


def render_manhattan_points(
    backend: PlotBackend,
    ax: Any,
    prepared_df: pd.DataFrame,
    chrom_order: Sequence[str],
    *,
    hover: Optional[HoverConfig] = None,
) -> None:
    """Scatter one prepared Manhattan panel, one chromosome colour at a time.

    When ``hover`` is supplied and the backend supports it, per-chromosome hover
    data is attached so interactive backends can show SNP tooltips.
    """
    for chrom in chrom_order:
        chrom_data = prepared_df[prepared_df["_chrom_str"] == chrom]
        if chrom_data.empty:
            continue
        hover_data = None
        if hover is not None and backend.supports_hover:
            hover_data = HoverDataBuilder(hover).build_dataframe(chrom_data)
        backend.scatter(
            ax,
            chrom_data["_cumulative_pos"],
            chrom_data["_neg_log_p"],
            colors=chrom_data["_color"].iloc[0],
            sizes=MANHATTAN_POINT_SIZE,
            marker="o",
            edgecolor=POINT_EDGE_COLOR,
            linewidth=MANHATTAN_EDGE_WIDTH,
            zorder=2,
            hover_data=hover_data,
        )

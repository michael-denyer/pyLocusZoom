"""Presentation policy more than one panel draws with."""

from typing import Any, Literal, Optional

import numpy as np

from ..backends.base import PlotBackend

REGIONAL_LINE_ALPHA = 0.65

# Manhattan and QQ styling constants
MANHATTAN_POINT_SIZE = 10
MANHATTAN_CATEGORICAL_POINT_SIZE = 30
QQ_POINT_SIZE = 10
POINT_EDGE_COLOR = "black"
MANHATTAN_EDGE_WIDTH = 0.1
QQ_EDGE_WIDTH = 0.02
QQ_CI_ALPHA = 0.5
SIGNIFICANCE_LINE_COLOR = "red"
SUGGESTIVE_LINE_COLOR = "blue"


def add_significance_line(
    backend: PlotBackend,
    ax: Any,
    threshold: Optional[float],
    *,
    axis: Literal["x", "y"] = "y",
    color: str = SIGNIFICANCE_LINE_COLOR,
    alpha: float = 1.0,
    linestyle: str = "--",
    linewidth: float = 1,
) -> None:
    """Draw the significance line at ``-log10(threshold)``.

    Args:
        backend: Plot backend instance.
        ax: Axes object from backend.
        threshold: P-value threshold (e.g., 5e-8). None to skip.
        axis: The axis the p-value is plotted on: ``"y"`` draws a horizontal
            line, ``"x"`` a vertical one.
        color: Line colour.
        alpha: Opacity of the line.
        linestyle: Matplotlib linestyle of the line.
        linewidth: Width of the line.
    """
    if threshold is None:
        return
    value = -np.log10(threshold)
    style = dict(
        color=color, linestyle=linestyle, linewidth=linewidth, alpha=alpha, zorder=1
    )
    if axis == "x":
        backend.axvline(ax, x=value, **style)
    else:
        backend.axhline(ax, y=value, **style)

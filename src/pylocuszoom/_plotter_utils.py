"""Shared utilities for plotter classes.

Internal module - not part of public API.
"""

from typing import Any, Literal, Optional, Union

import numpy as np

# Significance thresholds
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8
DEFAULT_EQTL_THRESHOLD = 1e-5


class _Unset:
    """Type of the ``UNSET`` sentinel."""

    def __repr__(self) -> str:
        return "UNSET"


# Marks a threshold argument the caller did not supply, so it falls back to the
# plotter's genomewide_threshold. None cannot carry that meaning: it already
# means "draw no significance line", and a caller must keep being able to say so.
UNSET = _Unset()

# A significance threshold as a caller may express it: a p-value, None for no
# line, or UNSET to inherit the plotter's own threshold.
ThresholdArg = Union[float, None, _Unset]


def resolve_threshold(
    supplied: ThresholdArg, plotter_default: Optional[float]
) -> Optional[float]:
    """Pick the threshold to draw at, honouring an explicit None.

    Args:
        supplied: The caller's threshold argument.
        plotter_default: The plotter's own ``genomewide_threshold``.

    Returns:
        The p-value to draw the line at, or None to draw no line.
    """
    return plotter_default if isinstance(supplied, _Unset) else supplied


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


def add_significance_line(
    backend: Any,
    ax: Any,
    threshold: Optional[float],
    *,
    axis: Literal["x", "y"] = "y",
    color: str = SIGNIFICANCE_LINE_COLOR,
    alpha: float = 1.0,
) -> None:
    """Draw the dashed significance line at ``-log10(threshold)``.

    Args:
        backend: Plot backend instance.
        ax: Axes object from backend.
        threshold: P-value threshold (e.g., 5e-8). None to skip.
        axis: The axis the p-value is plotted on: ``"y"`` draws a horizontal
            line, ``"x"`` a vertical one.
        color: Line colour.
        alpha: Opacity of the line.
    """
    if threshold is None:
        return
    value = -np.log10(threshold)
    style = dict(color=color, linestyle="--", linewidth=1, alpha=alpha, zorder=1)
    if axis == "x":
        backend.axvline(ax, x=value, **style)
    else:
        backend.axhline(ax, y=value, **style)

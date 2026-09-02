"""One QQ panel: the typed request and the function that draws it."""

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from ._plotter_utils import (
    POINT_EDGE_COLOR,
    QQ_CI_ALPHA,
    QQ_EDGE_WIDTH,
    QQ_POINT_SIZE,
    SIGNIFICANCE_LINE_COLOR,
)
from .backends.base import PlotBackend
from .colors import QQ_CI_COLOR, QQ_POINT_COLOR


def qq_title(lambda_gc: float, *, show_lambda: bool, compact: bool) -> str:
    """Return the title a QQ panel carries when the caller names none.

    Args:
        lambda_gc: Genomic inflation factor from ``prepare_qq_data``.
        show_lambda: Whether to name the inflation factor.
        compact: Whether the panel is one of a stack, which has room for a
            shorter title.

    Returns:
        The title.
    """
    if show_lambda:
        return f"λ = {lambda_gc:.3f}" if compact else f"QQ Plot (λ = {lambda_gc:.3f})"
    return "QQ" if compact else "QQ Plot"


@dataclass(frozen=True)
class QQPanelSpec:
    """One QQ panel's data and presentation policy.

    Attributes:
        qq_df: Frame from ``prepare_qq_data``.
        show_confidence_band: Whether to shade the 95% band.
        title: Panel title.
        title_fontsize: Panel title size.
        label_fontsize: Axis label size.
        x_label: X axis label, or None for none.
        y_label: Y axis label.
    """

    qq_df: pd.DataFrame
    show_confidence_band: bool
    title: str
    title_fontsize: int
    label_fontsize: int = 12
    x_label: Optional[str] = r"Expected $-\log_{10}(p)$"
    y_label: str = r"Observed $-\log_{10}(p)$"

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw this panel onto a backend axis."""
        render_qq_panel(backend, ax, self)


def render_qq_panel(backend: PlotBackend, ax: Any, spec: QQPanelSpec) -> None:
    """Draw one QQ panel onto a backend axis.

    Args:
        backend: Backend that owns the drawing primitives.
        ax: Axis to draw onto.
        spec: The panel's data and presentation policy.
    """
    qq_df = spec.qq_df
    if spec.show_confidence_band:
        backend.fill_between(
            ax,
            x=qq_df["_expected"],
            y1=qq_df["_ci_lower"],
            y2=qq_df["_ci_upper"],
            color=QQ_CI_COLOR,
            alpha=QQ_CI_ALPHA,
            zorder=1,
        )

    max_val = max(qq_df["_expected"].max(), qq_df["_observed"].max())
    backend.line(
        ax,
        x=pd.Series([0, max_val]),
        y=pd.Series([0, max_val]),
        color=SIGNIFICANCE_LINE_COLOR,
        linestyle="--",
        linewidth=1,
        zorder=2,
    )
    backend.scatter(
        ax,
        qq_df["_expected"],
        qq_df["_observed"],
        colors=QQ_POINT_COLOR,
        sizes=QQ_POINT_SIZE,
        marker="o",
        edgecolor=POINT_EDGE_COLOR,
        linewidth=QQ_EDGE_WIDTH,
        zorder=3,
    )
    backend.set_xlim(ax, 0, max_val * 1.05)
    backend.set_ylim(ax, 0, max_val * 1.05)
    if spec.x_label is not None:
        backend.set_xlabel(ax, spec.x_label, fontsize=spec.label_fontsize)
    backend.set_ylabel(ax, spec.y_label, fontsize=spec.label_fontsize)
    backend.set_title(ax, spec.title, fontsize=spec.title_fontsize)

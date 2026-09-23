"""One QQ panel: the typed request, which draws itself."""

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from ..backends.base import PlotBackend
from ..colors import QQ_CI_COLOR, QQ_POINT_COLOR
from ..config import GenomeWideStyle
from ._shared import (
    POINT_EDGE_COLOR,
    QQ_CI_ALPHA,
    QQ_EDGE_WIDTH,
    QQ_POINT_SIZE,
    SIGNIFICANCE_LINE_COLOR,
)
from .manhattan import styled


def qq_title(lambda_gc: float, *, show_lambda: bool, compact: bool) -> Optional[str]:
    """Return the title a QQ panel carries when the caller names none.

    Args:
        lambda_gc: Genomic inflation factor from ``prepare_qq_data``.
        show_lambda: Whether to name the inflation factor.
        compact: Whether the panel sits beside a Manhattan panel, where the
            axes already say it is a QQ plot and only the inflation factor
            is worth a title.

    Returns:
        The title, or None for none.
    """
    if show_lambda:
        return f"λ = {lambda_gc:.3f}" if compact else f"QQ Plot (λ = {lambda_gc:.3f})"
    return None if compact else "QQ Plot"


@dataclass(frozen=True)
class QQPanelSpec:
    """One QQ panel's data and presentation policy.

    Attributes:
        qq_df: Frame from ``prepare_qq_data``.
        show_confidence_band: Whether to shade the 95% band.
        title: Panel title, or None for none.
        title_fontsize: Panel title size.
        label_fontsize: Axis label size.
        x_label: X axis label, or None for none.
        y_label: Y axis label.
        style: Caller styling. A field it sets overrides the matching
            field above.
    """

    qq_df: pd.DataFrame
    show_confidence_band: bool
    title: Optional[str]
    title_fontsize: int
    label_fontsize: int = 12
    x_label: Optional[str] = r"Expected $-\log_{10}(p)$"
    y_label: str = r"Observed $-\log_{10}(p)$"
    style: GenomeWideStyle = GenomeWideStyle()

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw this panel onto a backend axis."""
        qq_df = self.qq_df
        style = self.style
        if self.show_confidence_band:
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
            linestyle=style.line_style,
            linewidth=style.line_width,
            zorder=2,
        )
        backend.scatter(
            ax,
            qq_df["_expected"],
            qq_df["_observed"],
            colors=QQ_POINT_COLOR,
            sizes=styled(style.point_size, QQ_POINT_SIZE),
            marker="o",
            edgecolor=POINT_EDGE_COLOR,
            linewidth=styled(style.point_edge_width, QQ_EDGE_WIDTH),
            zorder=3,
            alpha=style.point_alpha,
        )
        backend.set_xlim(ax, 0, max_val * 1.05)
        backend.set_ylim(ax, 0, max_val * 1.05)
        if style.tick_label_fontsize is not None:
            backend.set_tick_fontsize(ax, style.tick_label_fontsize)
        label_fontsize = styled(style.axis_label_fontsize, self.label_fontsize)
        if self.x_label is not None:
            backend.set_xlabel(ax, self.x_label, fontsize=label_fontsize)
        backend.set_ylabel(ax, self.y_label, fontsize=label_fontsize)
        if self.title:
            backend.set_title(
                ax,
                self.title,
                fontsize=styled(style.panel_title_fontsize, self.title_fontsize),
                fontweight=style.title_fontweight,
            )

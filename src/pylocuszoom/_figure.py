"""One figure plan for every plot family, and the one function that renders it.

A figure is an ordered list of panels laid out on a grid, plus the policy
that applies to the figure as a whole: its size, the label on the shared x
axis, spans highlighted across every panel, a title, and the layout
fractions.  Each plotter validates its input, prepares one panel value per
cell, and hands a :class:`FigurePlan` to :func:`render_figure`, which is the
only code above ``backends/`` that creates a figure or finalizes its layout.
"""

from dataclasses import dataclass
from typing import Any, Optional, Protocol, Sequence, Tuple

from .backends.base import PlotBackend


class Panel(Protocol):
    """A prepared panel that draws itself onto one backend axis."""

    def draw(self, backend: PlotBackend, ax: Any) -> None: ...


@dataclass(frozen=True)
class RegionHighlight:
    """A vertical span drawn across every panel, in x-axis coordinates."""

    start: float
    end: float
    color: str
    alpha: float


@dataclass(frozen=True)
class FigurePlan:
    """Everything one figure needs beyond what its panels draw.

    Attributes:
        panels: Panels in row-major cell order.
        figsize: Figure size as (width, height) in inches.
        height_ratios: Relative row heights, or None for equal rows.
        n_cols: Columns in the grid; one gives a vertical stack.
        width_ratios: Relative column widths for a grid, or None for equal.
        sharex: Whether stacked panels share their x axis.
        xlabel: Label for the bottom panel's x axis, or None for none.
        mb_xaxis: Whether every x axis reads in megabases.
        highlights: Spans highlighted across every panel.
        suptitle: Figure-level title, or None for none.
        top: Fraction of the figure height the panels extend to.
        hspace: Vertical space between panels as a fraction of panel height.
    """

    panels: Sequence[Panel]
    figsize: Tuple[float, float]
    height_ratios: Optional[Sequence[float]] = None
    n_cols: int = 1
    width_ratios: Optional[Sequence[float]] = None
    sharex: bool = True
    xlabel: Optional[str] = None
    mb_xaxis: bool = False
    highlights: Sequence[RegionHighlight] = ()
    suptitle: Optional[str] = None
    top: float = 0.95
    hspace: float = 0.08


def render_figure(backend: PlotBackend, plan: FigurePlan) -> Any:
    """Create the figure a plan describes, draw its panels, and finalize it.

    Raises:
        ValueError: If the plan has no panels.
    """
    n_panels = len(plan.panels)
    if n_panels == 0:
        raise ValueError("Figure plan must contain at least one panel")

    if plan.n_cols == 1:
        height_ratios = plan.height_ratios or [1.0] * n_panels
        fig, axes = backend.create_figure(
            height_ratios=list(height_ratios),
            figsize=plan.figsize,
            sharex=plan.sharex,
        )
    else:
        fig, axes = backend.create_figure_grid(
            n_rows=n_panels // plan.n_cols,
            n_cols=plan.n_cols,
            width_ratios=None if plan.width_ratios is None else list(plan.width_ratios),
            height_ratios=None
            if plan.height_ratios is None
            else list(plan.height_ratios),
            figsize=plan.figsize,
        )

    for ax, panel in zip(axes, plan.panels):
        panel.draw(backend, ax)

    # The label goes on before the tick format: plotly writes layout keys in
    # insertion order, so swapping these two changes every regional export.
    if plan.xlabel is not None:
        backend.set_xlabel(axes[-1], plan.xlabel)
    if plan.mb_xaxis:
        for ax in axes:
            backend.format_xaxis_mb(ax)
    for span in plan.highlights:
        backend.add_region_highlight(
            axes, span.start, span.end, color=span.color, alpha=span.alpha
        )
    if plan.suptitle:
        backend.set_suptitle(fig, plan.suptitle, fontsize=14)
    backend.finalize_layout(fig, top=plan.top, hspace=plan.hspace)
    return fig

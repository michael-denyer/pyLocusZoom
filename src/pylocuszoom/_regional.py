"""Shared regional association composition.

This module owns the figure-level policy shared by single and stacked regional
plots: figure creation, axis assignment, and the shared genomic x axis.  The
panels themselves and their drawing live in ``_regional_panels.py``.
"""

from typing import Any

from ._regional_panels import RegionalFigurePlan
from .backends.base import PlotBackend


def render_regional(backend: PlotBackend, plan: RegionalFigurePlan) -> Any:
    """Render every panel in one figure plan and finalize its layout."""
    if not plan.panels:
        raise ValueError("Regional figure plan must contain at least one panel")

    fig, axes = backend.create_figure(
        height_ratios=[panel.height for panel in plan.panels],
        figsize=plan.figsize,
        sharex=True,
    )
    for ax, panel in zip(axes, plan.panels):
        panel.draw(backend, ax)

    backend.set_xlabel(axes[-1], f"Chromosome {plan.chrom} (Mb)")
    for ax in axes:
        backend.format_xaxis_mb(ax)
    backend.finalize_layout(fig, hspace=plan.hspace)
    return fig

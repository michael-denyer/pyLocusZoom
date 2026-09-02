"""Shared regional association composition.

This module owns the figure-level policy shared by single and stacked regional
plots: figure creation, axis assignment, and dispatch to the drawing function
for each panel type.  The panels themselves and their drawing live in
``_regional_panels.py``.
"""

from functools import singledispatchmethod
from typing import Any

from ._regional_panels import (
    AssociationPanel,
    EqtlPanel,
    FinemappingPanel,
    GenePanel,
    HeatmapPanel,
    RegionalFigurePlan,
    RegionalPanel,
    draw_association,
    draw_eqtl,
    draw_finemapping,
    draw_genes,
    draw_heatmap,
)
from .backends.base import PlotBackend


class RegionalPlotComposer:
    """Compose shared association-panel policy for regional plot modes."""

    def __init__(
        self,
        backend: PlotBackend,
        genomewide_threshold: float,
    ):
        self._backend = backend
        self._genomewide_threshold = genomewide_threshold

    def render(self, plan: RegionalFigurePlan) -> Any:
        """Render every panel in one figure plan and finalize its layout."""
        if not plan.panels:
            raise ValueError("Regional figure plan must contain at least one panel")

        fig, axes = self._backend.create_figure(
            height_ratios=[panel.height for panel in plan.panels],
            figsize=plan.figsize,
            sharex=True,
        )
        for ax, panel in zip(axes, plan.panels):
            self.render_panel(panel, ax, fig, plan)

        self._backend.set_xlabel(axes[-1], f"Chromosome {plan.chrom} (Mb)")
        for ax in axes:
            self._backend.format_xaxis_mb(ax)
        self._backend.finalize_layout(fig, hspace=plan.hspace)
        return fig

    @singledispatchmethod
    def render_panel(
        self, panel: RegionalPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        """Render one panel onto its axes, dispatching on the panel type."""
        raise TypeError(f"No renderer for {type(panel).__name__}")

    @render_panel.register
    def _render_association(
        self, panel: AssociationPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        draw_association(self._backend, ax, panel, plan, self._genomewide_threshold)

    @render_panel.register
    def _render_heatmap(
        self, panel: HeatmapPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        draw_heatmap(self._backend, ax, panel, plan)

    @render_panel.register
    def _render_genes(
        self, panel: GenePanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        draw_genes(self._backend, ax, panel, plan)

    @render_panel.register
    def _render_finemapping(
        self, panel: FinemappingPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        draw_finemapping(self._backend, ax, panel)

    @render_panel.register
    def _render_eqtl(
        self, panel: EqtlPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        draw_eqtl(self._backend, ax, panel)

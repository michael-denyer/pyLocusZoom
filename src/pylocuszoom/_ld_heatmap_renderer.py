"""Semantic renderer for standalone LD heatmaps."""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import numpy as np

from .backends.base import PlotBackend, SupportsHeatmap
from .backends.composition import heatmap_highlight_rects, lower_triangle
from .colors import (
    LD_HEATMAP_COLORS,
    LEAD_SNP_HIGHLIGHT_COLOR,
    SECONDARY_HIGHLIGHT_COLOR,
)


@dataclass(frozen=True)
class LDHeatmapRequest:
    """One standalone LD heatmap, resolved by the plotter.

    ``lead_idx`` and ``highlight_indices`` index into ``snp_ids``.
    """

    data: np.ndarray
    snp_ids: List[str]
    lead_idx: Optional[int]
    highlight_indices: List[int]
    metric: str
    figsize: Tuple[float, float]
    title: Optional[str]
    show_colorbar: bool


def require_heatmap_backend(backend: PlotBackend) -> SupportsHeatmap:
    """Return the backend narrowed to ``SupportsHeatmap``.

    Raises:
        TypeError: If the backend does not implement ``SupportsHeatmap``. The
            heatmap is the whole figure here, so there is nothing to degrade to.
    """
    if not isinstance(backend, SupportsHeatmap):
        raise TypeError(
            f"{type(backend).__name__} does not support heatmaps. "
            "An LD heatmap needs a backend implementing SupportsHeatmap "
            "(add_heatmap, add_colorbar)."
        )
    return backend


def render_ld_heatmap(backend: SupportsHeatmap, req: LDHeatmapRequest) -> Any:
    """Draw the lower-triangle heatmap, its highlights, ticks, and title."""
    n_snps = len(req.snp_ids)
    fig, axes = backend.create_figure(
        height_ratios=[1.0], figsize=req.figsize, sharex=False
    )
    ax = axes[0]
    mappable = backend.add_heatmap(
        ax,
        data=lower_triangle(req.data),
        x_coords=list(range(n_snps)),
        y_coords=list(range(n_snps)),
        cmap_colors=LD_HEATMAP_COLORS,
        vmin=0.0,
        vmax=1.0,
    )
    if req.show_colorbar:
        backend.add_colorbar(ax, mappable, label="R²" if req.metric == "r2" else "D'")
    if req.lead_idx is not None:
        _highlight(backend, ax, req.lead_idx, n_snps, LEAD_SNP_HIGHLIGHT_COLOR)
    for idx in req.highlight_indices:
        _highlight(backend, ax, idx, n_snps, SECONDARY_HIGHLIGHT_COLOR)
    ticks = list(range(n_snps))
    backend.set_xticks(ax, ticks, req.snp_ids, rotation=90)
    backend.set_yticks(ax, ticks, req.snp_ids)
    if req.title:
        backend.set_title(ax, req.title)
    backend.finalize_layout(fig)
    return fig


def _highlight(
    backend: PlotBackend, ax: Any, idx: int, n_snps: int, color: str
) -> None:
    coords = list(range(n_snps))
    for x0, y0, width, height in heatmap_highlight_rects(idx, coords, coords):
        backend.add_rectangle(
            ax,
            (x0, y0),
            width,
            height,
            facecolor=None,
            edgecolor=color,
            linewidth=2,
            zorder=10,
        )

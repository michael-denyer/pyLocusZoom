"""The standalone LD heatmap panel, which draws itself."""

from dataclasses import dataclass
from typing import Any, List, Optional

import numpy as np

from .backends.base import PlotBackend
from .backends.composition import heatmap_highlight_rects, lower_triangle
from .colors import (
    LD_HEATMAP_COLORS,
    LEAD_SNP_HIGHLIGHT_COLOR,
    SECONDARY_HIGHLIGHT_COLOR,
)


@dataclass(frozen=True)
class LDHeatmapPanel:
    """One standalone LD heatmap, resolved by the plotter.

    ``lead_idx`` and ``highlight_indices`` index into ``snp_ids``.
    """

    data: np.ndarray
    snp_ids: List[str]
    lead_idx: Optional[int]
    highlight_indices: List[int]
    metric: str
    title: Optional[str]
    show_colorbar: bool

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the lower-triangle heatmap, its highlights, ticks, and title."""
        n_snps = len(self.snp_ids)
        mappable = backend.add_heatmap(
            ax,
            data=lower_triangle(self.data),
            x_coords=list(range(n_snps)),
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
        )
        if self.show_colorbar:
            backend.add_colorbar(
                ax, mappable, label="R²" if self.metric == "r2" else "D'"
            )
        if self.lead_idx is not None:
            _highlight(backend, ax, self.lead_idx, n_snps, LEAD_SNP_HIGHLIGHT_COLOR)
        for idx in self.highlight_indices:
            _highlight(backend, ax, idx, n_snps, SECONDARY_HIGHLIGHT_COLOR)
        ticks = list(range(n_snps))
        backend.set_xticks(ax, ticks, self.snp_ids, rotation=90)
        backend.set_yticks(ax, ticks, self.snp_ids)
        if self.title:
            backend.set_title(ax, self.title)


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

"""Semantic renderer for standalone LD heatmaps."""

from typing import Any, List, Optional, Tuple

import numpy as np

from .backends.base import PlotBackend, SupportsHeatmap
from .backends.composition import heatmap_highlight_rects, lower_triangle
from .colors import (
    LD_HEATMAP_COLORS,
    LEAD_SNP_HIGHLIGHT_COLOR,
    SECONDARY_HIGHLIGHT_COLOR,
)


class LDHeatmapRenderer:
    """Render a prepared LD matrix intent.

    Raises:
        TypeError: If the backend does not implement ``SupportsHeatmap``. The
            heatmap is the whole figure here, so there is nothing to degrade to.
    """

    def __init__(self, backend: PlotBackend):
        if not isinstance(backend, SupportsHeatmap):
            raise TypeError(
                f"{type(backend).__name__} does not support heatmaps. "
                "An LD heatmap needs a backend implementing SupportsHeatmap "
                "(add_heatmap, add_colorbar)."
            )
        self._backend = backend

    def render(
        self,
        data: np.ndarray,
        snp_ids: List[str],
        *,
        lead_idx: Optional[int],
        highlight_indices: List[int],
        metric: str,
        figsize: Tuple[float, float],
        title: Optional[str],
        show_colorbar: bool,
    ) -> Any:
        n_snps = len(snp_ids)
        fig, axes = self._backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=figsize, sharex=False
        )
        ax = axes[0]
        mappable = self._backend.add_heatmap(
            ax,
            data=lower_triangle(data),
            x_coords=list(range(n_snps)),
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
        )
        if show_colorbar:
            self._backend.add_colorbar(
                ax, mappable, label="R²" if metric == "r2" else "D'"
            )
        if lead_idx is not None:
            self._highlight(ax, lead_idx, n_snps, LEAD_SNP_HIGHLIGHT_COLOR)
        for idx in highlight_indices:
            self._highlight(ax, idx, n_snps, SECONDARY_HIGHLIGHT_COLOR)
        ticks = list(range(n_snps))
        self._backend.set_xticks(ax, ticks, snp_ids, rotation=90)
        self._backend.set_yticks(ax, ticks, snp_ids)
        if title:
            self._backend.set_title(ax, title)
        self._backend.finalize_layout(fig)
        return fig

    def _highlight(self, ax: Any, idx: int, n_snps: int, color: str) -> None:
        coords = list(range(n_snps))
        for x0, y0, width, height in heatmap_highlight_rects(idx, coords, coords):
            self._backend.add_rectangle(
                ax,
                (x0, y0),
                width,
                height,
                facecolor=None,
                edgecolor=color,
                linewidth=2,
                zorder=10,
            )

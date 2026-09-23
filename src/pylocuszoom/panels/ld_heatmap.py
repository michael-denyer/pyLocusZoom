"""The standalone LD heatmap panel, which draws itself."""

from dataclasses import dataclass
from typing import Any, List, Optional

import numpy as np

from ..backends.base import PlotBackend
from ..backends.composition import draw_ld_heatmap
from ..colors import LEAD_SNP_HIGHLIGHT_COLOR, SECONDARY_HIGHLIGHT_COLOR


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
        ticks = list(range(len(self.snp_ids)))
        lead = [] if self.lead_idx is None else [self.lead_idx]
        draw_ld_heatmap(
            backend,
            ax,
            self.data,
            ticks,
            metric=self.metric,
            show_colorbar=self.show_colorbar,
            outlines=[(idx, LEAD_SNP_HIGHLIGHT_COLOR) for idx in lead]
            + [(idx, SECONDARY_HIGHLIGHT_COLOR) for idx in self.highlight_indices],
        )
        backend.set_xticks(ax, ticks, self.snp_ids, rotation=90)
        backend.set_yticks(ax, ticks, self.snp_ids)
        if self.title:
            backend.set_title(ax, self.title)

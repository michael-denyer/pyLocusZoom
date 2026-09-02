"""The regional LD heatmap panel, drawn under an association panel."""

from dataclasses import dataclass
from typing import Any, List, Optional

import pandas as pd

from ..backends.base import PlotBackend
from ..backends.composition import heatmap_highlight_rects, lower_triangle
from ..colors import LD_HEATMAP_COLORS, LEAD_SNP_HIGHLIGHT_COLOR
from ..config import RegionConfig
from ..logging import logger
from .association import AssociationPanel


@dataclass(frozen=True)
class HeatmapPanel:
    """Prepared regional LD heatmap panel."""

    matrix: pd.DataFrame
    region: RegionConfig
    height: float
    x_positions: List[int]
    snp_ids: List[str]
    metric: str
    lead_snp_id: Optional[str]

    @classmethod
    def from_matrix(
        cls,
        ld_matrix: pd.DataFrame,
        snp_ids: List[str],
        *,
        source: AssociationPanel,
        region: RegionConfig,
        height: float,
        metric: str,
    ) -> "HeatmapPanel":
        """Map heatmap SNP ids to positions through the source panel's frame.

        Raises:
            ValueError: If the source frame has no SNP id column, or no
                heatmap SNP falls inside the region.
        """
        df = source.data
        rs_col, pos_col = source.columns.rs_col, source.columns.pos_col
        if rs_col not in df.columns:
            raise ValueError(
                f"Cannot map heatmap to genomic coords: column '{rs_col}' not in GWAS data"
            )

        snp_to_pos = dict(zip(df[rs_col], df[pos_col]))
        kept = [
            (i, snp_id, int(snp_to_pos[snp_id]))
            for i, snp_id in enumerate(snp_ids)
            if snp_id in snp_to_pos and region.start <= snp_to_pos[snp_id] <= region.end
        ]
        if not kept:
            raise ValueError(
                "No SNPs from LD heatmap overlap with region - heatmap not rendered"
            )
        indices, kept_ids, x_positions = (list(column) for column in zip(*kept))

        lead_snp_id = None
        if source.lead_pos is not None:
            lead_row = df[df[pos_col] == source.lead_pos]
            if not lead_row.empty:
                lead_snp_id = lead_row[rs_col].iloc[0]
        return cls(
            matrix=ld_matrix.iloc[indices, indices].copy(),
            region=region,
            height=height,
            x_positions=x_positions,
            snp_ids=kept_ids,
            metric=metric,
            lead_snp_id=lead_snp_id,
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the lower-triangle LD heatmap and its lead-SNP crosshair."""
        n_snps = len(self.snp_ids)
        if n_snps < 2:
            logger.debug("Skipping heatmap: fewer than 2 SNPs after filtering")
            return
        mappable = backend.add_heatmap(
            ax,
            data=lower_triangle(self.matrix.values),
            x_coords=self.x_positions,
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
        )
        backend.add_colorbar(ax, mappable, label="R²" if self.metric == "r2" else "D'")
        if self.lead_snp_id is not None and self.lead_snp_id in self.snp_ids:
            rects = heatmap_highlight_rects(
                self.snp_ids.index(self.lead_snp_id),
                self.x_positions,
                list(range(n_snps)),
            )
            for x0, y0, width, height in rects:
                backend.add_rectangle(
                    ax,
                    (x0, y0),
                    width,
                    height,
                    facecolor=None,
                    edgecolor=LEAD_SNP_HIGHLIGHT_COLOR,
                    linewidth=2,
                    zorder=10,
                )
        backend.set_xlim(ax, self.region.start, self.region.end)
        backend.hide_yaxis(ax)

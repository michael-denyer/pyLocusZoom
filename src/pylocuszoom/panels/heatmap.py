"""The regional LD heatmap panel, drawn under an association panel."""

from dataclasses import dataclass
from typing import Any, List, Optional

import pandas as pd

from ..backends.base import PlotBackend
from ..backends.composition import draw_ld_heatmap
from ..colors import LEAD_SNP_HIGHLIGHT_COLOR
from ..config import RegionConfig
from ..exceptions import ValidationError
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
            ValidationError: If the source frame has no SNP id column, or no
                heatmap SNP falls inside the region.
        """
        df = source.data
        rs_col, pos_col = source.hover.snp_col, source.columns.pos_col
        if rs_col is None:
            raise ValidationError(
                "Cannot map heatmap to genomic coords: column "
                f"'{source.columns.rs_col}' not in GWAS data"
            )

        snp_to_pos = dict(zip(df[rs_col], df[pos_col]))
        kept = [
            (i, snp_id, int(snp_to_pos[snp_id]))
            for i, snp_id in enumerate(snp_ids)
            if snp_id in snp_to_pos and region.start <= snp_to_pos[snp_id] <= region.end
        ]
        if not kept:
            raise ValidationError(
                "No SNPs from LD heatmap overlap with region - heatmap not rendered"
            )
        kept.sort(key=lambda record: record[2])
        indices, kept_ids, x_positions = (list(column) for column in zip(*kept))
        if len(set(x_positions)) != len(x_positions):
            raise ValidationError(
                "Regional heatmap SNPs must have distinct genomic positions"
            )

        lead_snp_id = (
            df.at[source.lead_index, rs_col] if source.lead_index is not None else None
        )
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
        draw_ld_heatmap(
            backend,
            ax,
            self.matrix.values,
            self.x_positions,
            metric=self.metric,
            show_colorbar=True,
            outlines=(
                [(self.snp_ids.index(self.lead_snp_id), LEAD_SNP_HIGHLIGHT_COLOR)]
                if self.lead_snp_id in self.snp_ids
                else []
            ),
        )
        backend.set_xlim(ax, self.region.start, self.region.end)
        backend.hide_yaxis(ax)

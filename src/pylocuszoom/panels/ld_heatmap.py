"""The standalone LD heatmap panel, which draws itself."""

from dataclasses import dataclass
from typing import Any, List, Optional, Union, get_args

import numpy as np
import pandas as pd

from .._ld_matrix import prepare_ld_matrix
from ..backends.base import PlotBackend
from ..backends.composition import draw_ld_heatmap
from ..colors import LEAD_SNP_HIGHLIGHT_COLOR, SECONDARY_HIGHLIGHT_COLOR
from ..config import LDMetric
from ..exceptions import ValidationError


@dataclass(frozen=True)
class LDHeatmapPanel:
    """One standalone LD heatmap, resolved when it is built.

    ``lead_idx`` and ``highlight_indices`` index into ``snp_ids``.
    """

    data: np.ndarray
    snp_ids: List[str]
    lead_idx: Optional[int]
    highlight_indices: List[int]
    metric: str
    title: Optional[str]
    show_colorbar: bool

    @classmethod
    def from_matrix(
        cls,
        ld_matrix: Union[pd.DataFrame, np.ndarray],
        snp_ids: Optional[List[str]],
        *,
        lead_snp: Optional[str],
        highlight_snps: Optional[List[str]],
        metric: str,
        title: Optional[str],
        show_colorbar: bool,
    ) -> "LDHeatmapPanel":
        """Validate a square LD matrix and resolve its SNP ids and highlights.

        ``snp_ids`` defaults to a DataFrame's index, or to ``"0" .. "n-1"``
        for an array.

        Raises:
            ValidationError: If the matrix is not square, ``snp_ids`` does not
                match it, ``lead_snp`` or a highlight SNP is not in
                ``snp_ids``, or ``metric`` is not ``"r2"`` or ``"dprime"``.
        """
        if metric not in get_args(LDMetric):
            raise ValidationError(f"metric must be 'r2' or 'dprime', got {metric!r}")

        data, snp_ids = prepare_ld_matrix(ld_matrix, snp_ids)

        if lead_snp is not None and lead_snp not in snp_ids:
            raise ValidationError(f"lead_snp '{lead_snp}' not found in snp_ids")
        for snp in highlight_snps or ():
            if snp not in snp_ids:
                raise ValidationError(f"highlight_snp '{snp}' not found in snp_ids")

        return cls(
            data=data,
            snp_ids=snp_ids,
            lead_idx=None if lead_snp is None else snp_ids.index(lead_snp),
            highlight_indices=[snp_ids.index(snp) for snp in highlight_snps or ()],
            metric=metric,
            title=title,
            show_colorbar=show_colorbar,
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the lower-triangle heatmap, its highlights, ticks, and title."""
        ticks = list(range(len(self.snp_ids)))
        outlines = [(idx, SECONDARY_HIGHLIGHT_COLOR) for idx in self.highlight_indices]
        if self.lead_idx is not None:
            outlines.insert(0, (self.lead_idx, LEAD_SNP_HIGHLIGHT_COLOR))
        draw_ld_heatmap(
            backend,
            ax,
            self.data,
            ticks,
            metric=self.metric,
            show_colorbar=self.show_colorbar,
            outlines=outlines,
        )
        backend.set_xticks(ax, ticks, self.snp_ids, rotation=90)
        backend.set_yticks(ax, ticks, self.snp_ids)
        if self.title:
            backend.set_title(ax, self.title)

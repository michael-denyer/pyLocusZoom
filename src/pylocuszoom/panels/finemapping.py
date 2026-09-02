"""The fine-mapping panel: the PIP line and its credible-set points."""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import pandas as pd

from ..backends.base import PlotBackend
from ..backends.composition import finemapping_legend_entries
from ..backends.hover import HoverConfig, HoverDataBuilder
from ..colors import NO_DATA_COLOR, PIP_LINE_COLOR, get_credible_set_color
from ..config import RegionConfig
from ..finemapping import get_credible_sets, prepare_finemapping_for_plotting

PIP_SCATTER_THRESHOLD = 0.01


@dataclass(frozen=True)
class FinemappingPanel:
    """Prepared fine-mapping panel.

    ``cs_col`` is a column of ``data`` whenever it is not None, and
    ``credible_sets`` are the set ids it holds.
    """

    data: pd.DataFrame
    height: float
    cs_col: Optional[str]
    credible_sets: List[int]
    hover: HoverConfig

    @classmethod
    def from_frame(
        cls, df: pd.DataFrame, region: RegionConfig, cs_col: Optional[str]
    ) -> "FinemappingPanel":
        """Validate, region-filter, and sort raw fine-mapping results."""
        data = prepare_finemapping_for_plotting(
            df,
            pos_col="pos",
            pip_col="pip",
            chrom=region.chrom,
            start=region.start,
            end=region.end,
        )
        resolved = cs_col if cs_col and cs_col in data.columns else None
        extra_cols = {"pip": "PIP"}
        if resolved:
            extra_cols[resolved] = "Credible Set"
        return cls(
            data=data,
            height=1.5,
            cs_col=resolved,
            credible_sets=get_credible_sets(data, resolved) if resolved else [],
            hover=HoverConfig(pos_col="pos", extra_cols=extra_cols),
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the PIP line, credible-set points, and their legend."""
        data = self.data
        if not data.empty:
            backend.line(
                ax,
                data["pos"],
                data["pip"],
                color=PIP_LINE_COLOR,
                linewidth=1.5,
                alpha=0.8,
                zorder=1,
            )
            hover_builder = HoverDataBuilder(self.hover)
            for subset, color, size, linewidth, zorder in _finemapping_groups(self):
                if subset.empty:
                    continue
                backend.scatter(
                    ax,
                    subset["pos"],
                    subset["pip"],
                    colors=color,
                    sizes=size,
                    marker="o",
                    edgecolor="black",
                    linewidth=linewidth,
                    zorder=zorder,
                    hover_data=hover_builder.build_dataframe(subset),
                )
            if self.credible_sets:
                backend.add_legend(
                    ax,
                    finemapping_legend_entries(self.credible_sets),
                    loc="upper right",
                    title="Credible sets",
                )
        backend.set_ylabel(ax, "PIP")
        backend.set_ylim(ax, -0.05, 1.05)


def _finemapping_groups(
    panel: FinemappingPanel,
) -> List[Tuple[pd.DataFrame, str, int, float, int]]:
    """Split the PIP points into scatter groups, each with its own styling."""
    data = panel.data
    if not panel.credible_sets:
        above = data[data["pip"] >= PIP_SCATTER_THRESHOLD]
        return [(above, PIP_LINE_COLOR, 50, 0.5, 3)]

    cs_values = data[panel.cs_col]
    groups = [
        (data[cs_values == cs_id], get_credible_set_color(cs_id), 50, 0.5, 3)
        for cs_id in panel.credible_sets
    ]
    unassigned = data[cs_values.isna() | (cs_values == 0)]
    groups.append(
        (
            unassigned[unassigned["pip"] >= PIP_SCATTER_THRESHOLD],
            NO_DATA_COLOR,
            30,
            0.3,
            2,
        )
    )
    return groups

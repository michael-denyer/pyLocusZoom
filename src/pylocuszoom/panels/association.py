"""The association panel: the scatter every regional figure is built around."""

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from .._plotter_utils import add_significance_line
from ..backends.base import (
    PlotBackend,
    SupportsSNPLabels,
)
from ..backends.composition import (
    LD_LEGEND_TITLE,
    ld_legend_entries,
    render_recombination_overlay,
)
from ..backends.hover import HoverConfig, HoverDataBuilder
from ..colors import (
    LEAD_SNP_COLOR,
    NO_DATA_COLOR,
    get_ld_bin,
    get_ld_color_palette,
)
from ..config import ColumnConfig, DisplayConfig, RegionConfig
from ._shared import REGIONAL_LINE_ALPHA


def hover_for_association(
    data: pd.DataFrame, columns: ColumnConfig, ld_col: Optional[str]
) -> HoverConfig:
    """Resolve the association hover contract against a prepared frame."""
    return HoverConfig(
        snp_col=columns.rs_col if columns.rs_col in data.columns else None,
        pos_col=columns.pos_col if columns.pos_col in data.columns else None,
        p_col=columns.p_col if columns.p_col in data.columns else None,
        ld_col=ld_col,
    )


@dataclass(frozen=True)
class AssociationPanel:
    """Prepared association panel and its presentation policy.

    ``data`` already carries ``neglog10p`` and any merged LD column.
    ``ld_col`` is a column of ``data`` whenever it is not None.
    """

    data: pd.DataFrame
    region: RegionConfig
    height: float
    columns: ColumnConfig
    display: DisplayConfig
    genomewide_threshold: Optional[float]
    ld_col: Optional[str]
    lead_pos: Optional[int]
    recomb_df: Optional[pd.DataFrame]
    hover: HoverConfig
    panel_label: Optional[str] = None
    add_ld_legend: bool = False

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the association scatter with its axes, overlay, and legends."""
        df = self.data
        columns = self.columns
        start, end = self.region.start, self.region.end
        _draw_association_points(backend, ax, self)
        add_significance_line(
            backend, ax, self.genomewide_threshold, alpha=REGIONAL_LINE_ALPHA
        )
        backend.set_ylabel(ax, r"$-\log_{10}$ P")
        y_max = df["neglog10p"].max()
        if pd.notna(y_max) and y_max > 0:
            backend.set_ylim(ax, 0, y_max * 1.15)
        backend.set_xlim(ax, start, end)

        if (
            self.display.snp_labels
            and columns.rs_col in df.columns
            and self.display.label_top_n > 0
            and not df.empty
            and isinstance(backend, SupportsSNPLabels)
        ):
            backend.add_snp_labels(
                ax,
                df,
                pos_col=columns.pos_col,
                neglog10p_col="neglog10p",
                rs_col=columns.rs_col,
                label_top_n=self.display.label_top_n,
                adjust=True,
                lead_pos=self.lead_pos,
                region_span=end - start,
            )

        recomb_df = self.recomb_df
        if recomb_df is not None and not recomb_df.empty:
            render_recombination_overlay(backend, ax, recomb_df, start, end)

        if self.panel_label:
            backend.add_panel_label(ax, self.panel_label)
        if self.add_ld_legend and self.ld_col is not None:
            backend.add_legend(
                ax, ld_legend_entries(), loc="upper right", title=LD_LEGEND_TITLE
            )


def _draw_association_points(
    backend: PlotBackend, ax: Any, panel: AssociationPanel
) -> None:
    """Draw association points, including LD and lead-SNP styling."""
    df = panel.data
    pos_col = panel.columns.pos_col
    ld_col = panel.ld_col
    hover_builder = HoverDataBuilder(panel.hover)

    if ld_col is not None:
        df = df.copy()
        df["ld_bin"] = df[ld_col].apply(get_ld_bin)
        df = df.sort_values(ld_col, ascending=True, na_position="first")
        palette = get_ld_color_palette()
        for bin_label in df["ld_bin"].unique():
            bin_data = df[df["ld_bin"] == bin_label]
            backend.scatter(
                ax,
                bin_data[pos_col],
                bin_data["neglog10p"],
                colors=palette[bin_label],
                sizes=60,
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
                hover_data=hover_builder.build_dataframe(bin_data),
            )
    else:
        backend.scatter(
            ax,
            df[pos_col],
            df["neglog10p"],
            colors=NO_DATA_COLOR,
            sizes=60,
            edgecolor="black",
            linewidth=0.5,
            zorder=2,
            hover_data=hover_builder.build_dataframe(df),
        )

    if panel.lead_pos is not None:
        lead_snp = df[df[pos_col] == panel.lead_pos]
        if not lead_snp.empty:
            backend.scatter(
                ax,
                lead_snp[pos_col],
                lead_snp["neglog10p"],
                colors=LEAD_SNP_COLOR,
                sizes=120,
                marker="D",
                edgecolor="black",
                linewidth=1.5,
                zorder=10,
                hover_data=hover_builder.build_dataframe(lead_snp),
            )

"""The association panel: the scatter every regional figure is built around."""

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from .._data import prepare_pvalue_data
from .._label_data import select_label_candidates
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
from ..config import ColumnConfig, DisplayConfig, LDConfig, RegionConfig
from ..exceptions import ValidationError
from ..logging import logger
from ..schemas import Canonical, gwas_plot_spec
from ..utils import filter_by_region
from ..validation import check, resolve_column
from ._shared import REGIONAL_LINE_ALPHA, add_significance_line


@dataclass(frozen=True)
class AssociationInput:
    """One region-selected association frame with its resolved column roles.

    ``rs_col`` and ``ld_col`` are columns of ``data`` whenever they are not
    None. LD enrichment replaces ``data`` and ``ld_col`` together.
    """

    data: pd.DataFrame
    columns: ColumnConfig
    rs_col: Optional[str]
    ld_col: Optional[str]
    ld_reference_file: Optional[str]
    lead_index: Optional[int]
    label: Optional[str] = None

    @classmethod
    def prepare(
        cls,
        frame: pd.DataFrame,
        region: RegionConfig,
        columns: ColumnConfig,
        ld: LDConfig,
        label: Optional[str] = None,
    ) -> "AssociationInput":
        """Validate a caller's frame, select the region and resolve the lead.

        Raises:
            ValidationError: If a named column is missing, or LD from a
                reference fileset has no SNP id column to look up.
        """
        check(frame, gwas_plot_spec(columns.pos_col, columns.p_col))
        resolve_column(frame, ld.ld_col, parameter="ld_col")
        rs_col = resolve_column(
            frame, columns.rs_col, parameter="rs_col", optional_default=Canonical.RS
        )
        if ld.ld_reference_file is not None and rs_col is None:
            raise ValidationError(
                "ld_reference_file needs SNP ids to compute LD, and column "
                f"'{columns.rs_col}' is not in the GWAS data. Add it, or name "
                "the SNP id column with ColumnConfig(rs_col=...)."
            )
        selected = filter_by_region(
            frame,
            region=(region.chrom, region.start, region.end),
            chrom_col=columns.chrom_col,
            pos_col=columns.pos_col,
        )
        data = prepare_pvalue_data(selected, columns.p_col, "regional")
        data = data.reset_index(drop=True)
        candidates = (
            data if ld.lead_pos is None else data[data[columns.pos_col] == ld.lead_pos]
        )
        lead_index = (
            int(candidates["neglog10p"].idxmax()) if not candidates.empty else None
        )
        if ld.lead_pos is not None and lead_index is None:
            logger.warning(
                "Lead SNP at position {} not found in region; LD coloring will be skipped",
                ld.lead_pos,
            )
        return cls(
            data, columns, rs_col, ld.ld_col, ld.ld_reference_file, lead_index, label
        )


@dataclass(frozen=True)
class AssociationPanel:
    """Prepared association panel and its presentation policy.

    ``data`` already carries ``neglog10p`` and any merged LD column.
    ``ld_col`` and ``hover.snp_col`` are columns of ``data`` whenever they
    are not None.
    """

    data: pd.DataFrame
    region: RegionConfig
    height: float
    columns: ColumnConfig
    display: DisplayConfig
    genomewide_threshold: Optional[float]
    ld_col: Optional[str]
    lead_index: Optional[int]
    recomb_df: Optional[pd.DataFrame]
    hover: HoverConfig
    panel_label: Optional[str] = None
    add_ld_legend: bool = False

    @classmethod
    def from_input(
        cls,
        request: AssociationInput,
        *,
        region: RegionConfig,
        display: DisplayConfig,
        threshold: Optional[float],
        height: float,
        recomb_df: Optional[pd.DataFrame],
        is_top: bool,
    ) -> "AssociationPanel":
        """Build a panel from a resolved input; the top panel carries the LD legend.

        Args:
            request: The region-selected frame and its column roles.
            region: The figure's region.
            display: Display options with the per-figure defaults applied.
            threshold: P-value for the significance line, or None.
            height: Height-ratio units of the panel.
            recomb_df: Recombination rates to overlay, or None.
            is_top: Whether this is the first association panel.
        """
        return cls(
            data=request.data,
            region=region,
            height=height,
            columns=request.columns,
            display=display,
            genomewide_threshold=threshold,
            ld_col=request.ld_col,
            lead_index=request.lead_index,
            recomb_df=recomb_df,
            hover=HoverConfig(
                snp_col=request.rs_col,
                pos_col=request.columns.pos_col,
                p_col=request.columns.p_col,
                ld_col=request.ld_col,
            ),
            panel_label=request.label,
            add_ld_legend=is_top,
        )

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

        snp_col = self.hover.snp_col
        if (
            self.display.snp_labels
            and snp_col is not None
            and self.display.label_top_n > 0
            and not df.empty
            and isinstance(backend, SupportsSNPLabels)
        ):
            backend.add_snp_labels(
                ax,
                select_label_candidates(
                    df,
                    pos_col=columns.pos_col,
                    lead_index=self.lead_index,
                    region_span=end - start,
                ),
                pos_col=columns.pos_col,
                neglog10p_col="neglog10p",
                rs_col=snp_col,
                label_top_n=self.display.label_top_n,
            )

        recomb_df = self.recomb_df
        if recomb_df is not None and not recomb_df.empty:
            render_recombination_overlay(backend, ax, recomb_df, start, end)

        if self.panel_label:
            backend.add_panel_label(ax, self.panel_label)
        if self.add_ld_legend and self.ld_col is not None:
            backend.add_legend(ax, ld_legend_entries(), title=LD_LEGEND_TITLE)


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
                hover_data=hover_builder.build(bin_data),
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
            hover_data=hover_builder.build(df),
        )

    if panel.lead_index is not None:
        lead_snp = df.loc[[panel.lead_index]]
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
            hover_data=hover_builder.build(lead_snp),
        )

"""Semantic renderer for Miami plots."""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._manhattan_panel import (
    chromosome_ticks,
    manhattan_spec,
    render_manhattan_panel,
    shared_manhattan_limits,
)
from .backends.base import PlotBackend, SupportsRegionHighlight
from .backends.hover import HoverConfig


class MiamiRenderer:
    """Render a prepared mirrored Manhattan figure."""

    def __init__(self, backend: PlotBackend):
        self._backend = backend

    def render(
        self,
        top_df: pd.DataFrame,
        bottom_df: pd.DataFrame,
        *,
        pos_col: str,
        p_col: str,
        rs_col: Optional[str],
        top_threshold: Optional[float],
        bottom_threshold: Optional[float],
        top_label: Optional[str],
        bottom_label: Optional[str],
        top_snp_annotations: Optional[List[str]],
        bottom_snp_annotations: Optional[List[str]],
        highlight_regions: Optional[List[Tuple[str, int, int]]],
        highlight_color: str,
        highlight_alpha: float,
        figsize: Tuple[float, float],
        title: Optional[str],
    ) -> Any:
        fig, axes = self._backend.create_figure(
            n_panels=2,
            height_ratios=[1.0, 1.0],
            figsize=figsize,
            sharex=True,
        )
        top_ax, bottom_ax = axes
        hover = (
            HoverConfig(snp_col=rs_col, pos_col=pos_col, p_col=p_col)
            if rs_col is not None
            else None
        )
        x_limits = shared_manhattan_limits([top_df, bottom_df])
        ticks = chromosome_ticks(
            top_df.attrs["chrom_order"], top_df.attrs["chrom_centers"]
        )
        render_manhattan_panel(
            self._backend,
            top_ax,
            manhattan_spec(
                top_df,
                x_limits=x_limits,
                ticks=ticks,
                significance_threshold=top_threshold,
                panel_label=top_label,
                hover=hover,
            ),
        )
        render_manhattan_panel(
            self._backend,
            bottom_ax,
            manhattan_spec(
                bottom_df,
                x_limits=x_limits,
                ticks=ticks,
                significance_threshold=bottom_threshold,
                x_label="Chromosome",
                panel_label=bottom_label,
                panel_label_y_frac=0.05,
                invert_y=True,
                hover=hover,
            ),
        )
        if top_snp_annotations and rs_col:
            self._add_snp_annotations(top_ax, top_df, rs_col, top_snp_annotations)
        if bottom_snp_annotations and rs_col:
            self._add_snp_annotations(
                bottom_ax, bottom_df, rs_col, bottom_snp_annotations
            )
        if highlight_regions:
            offsets = self._chrom_offsets(top_df, pos_col)
            for chrom, start, end in highlight_regions:
                self._add_region_highlight(
                    fig,
                    top_ax,
                    bottom_ax,
                    str(chrom),
                    start,
                    end,
                    offsets,
                    highlight_color,
                    highlight_alpha,
                )
        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.92, hspace=0.05)
        else:
            self._backend.finalize_layout(fig, hspace=0.05)
        return fig

    def _add_snp_annotations(
        self, ax: Any, prepared_df: pd.DataFrame, rs_col: str, snp_ids: List[str]
    ) -> None:
        for _, row in prepared_df[prepared_df[rs_col].isin(snp_ids)].iterrows():
            self._backend.add_text(
                ax,
                x=row["_cumulative_pos"],
                y=row["_neg_log_p"],
                text=str(row[rs_col]),
                fontsize=8,
                ha="center",
                va="bottom",
            )

    @staticmethod
    def _chrom_offsets(prepared_df: pd.DataFrame, pos_col: str) -> dict[str, float]:
        offsets = {}
        for chrom in prepared_df.attrs.get("chrom_order", []):
            chrom_data = prepared_df[prepared_df["_chrom_str"] == str(chrom)]
            if not chrom_data.empty:
                first_row = chrom_data.iloc[0]
                offsets[str(chrom)] = first_row["_cumulative_pos"] - first_row[pos_col]
        return offsets

    def _add_region_highlight(
        self,
        fig: Any,
        top_ax: Any,
        bottom_ax: Any,
        chrom: str,
        start: int,
        end: int,
        offsets: dict[str, float],
        color: str,
        alpha: float,
    ) -> None:
        if chrom not in offsets:
            return
        if not isinstance(self._backend, SupportsRegionHighlight):
            return
        self._backend.add_region_highlight(
            fig,
            [top_ax, bottom_ax],
            offsets[chrom] + start,
            offsets[chrom] + end,
            color=color,
            alpha=alpha,
        )

"""Semantic renderer for Miami plots."""

from typing import Any, List, Optional, Sequence, Tuple

import pandas as pd

from ._plotter_utils import (
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
    add_significance_line,
)
from .backends.base import PlotBackend, SupportsRegionHighlight
from .backends.hover import HoverConfig, HoverDataBuilder


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
        chrom_order = top_df.attrs["chrom_order"]
        chrom_centers = top_df.attrs["chrom_centers"]
        for ax, prepared_df, threshold in (
            (top_ax, top_df, top_threshold),
            (bottom_ax, bottom_df, bottom_threshold),
        ):
            self._render_points(ax, prepared_df, chrom_order, pos_col, p_col, rs_col)
            add_significance_line(self._backend, ax, threshold)

        x_min = min(top_df["_cumulative_pos"].min(), bottom_df["_cumulative_pos"].min())
        x_max = max(top_df["_cumulative_pos"].max(), bottom_df["_cumulative_pos"].max())
        x_padding = (x_max - x_min) * 0.01
        self._backend.set_xlim(top_ax, x_min - x_padding, x_max + x_padding)
        self._backend.set_xlim(bottom_ax, x_min - x_padding, x_max + x_padding)
        top_y_max = self._safe_ymax(top_df["_neg_log_p"].max())
        bottom_y_max = self._safe_ymax(bottom_df["_neg_log_p"].max())
        self._backend.set_ylim(top_ax, 0, top_y_max)
        self._backend.set_ylim(bottom_ax, bottom_y_max, 0)

        valid_chroms = [c for c in chrom_order if c in chrom_centers]
        positions = [chrom_centers[c] for c in valid_chroms]
        labels = [str(c) for c in valid_chroms]
        for ax in axes:
            self._backend.set_xticks(ax, positions, labels, fontsize=8)
        self._backend.set_ylabel(top_ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_ylabel(bottom_ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_xlabel(bottom_ax, "Chromosome", fontsize=12)
        if top_label:
            self._backend.add_panel_label(top_ax, top_label, y_frac=0.95)
        if bottom_label:
            self._backend.add_panel_label(bottom_ax, bottom_label, y_frac=0.05)
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
        self._backend.hide_spines(top_ax, ["top", "right"])
        self._backend.hide_spines(bottom_ax, ["top", "right"])
        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.92, hspace=0.05)
        else:
            self._backend.finalize_layout(fig, hspace=0.05)
        return fig

    @staticmethod
    def _safe_ymax(value: float) -> float:
        return max(value * 1.1, 1.0) if pd.notna(value) else 1.0

    def _render_points(
        self,
        ax: Any,
        prepared_df: pd.DataFrame,
        chrom_order: Sequence[str],
        pos_col: str,
        p_col: str,
        rs_col: Optional[str],
    ) -> None:
        for chrom in chrom_order:
            chrom_data = prepared_df[prepared_df["_chrom_str"] == chrom]
            if chrom_data.empty:
                continue
            hover_data = None
            if self._backend.supports_hover and rs_col is not None:
                hover_data = HoverDataBuilder(
                    HoverConfig(snp_col=rs_col, pos_col=pos_col, p_col=p_col)
                ).build_dataframe(chrom_data)
            self._backend.scatter(
                ax,
                chrom_data["_cumulative_pos"],
                chrom_data["_neg_log_p"],
                colors=chrom_data["_color"].iloc[0],
                sizes=MANHATTAN_POINT_SIZE,
                marker="o",
                edgecolor=POINT_EDGE_COLOR,
                linewidth=MANHATTAN_EDGE_WIDTH,
                zorder=2,
                hover_data=hover_data,
            )

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

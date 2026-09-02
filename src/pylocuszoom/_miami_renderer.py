"""Semantic renderer for Miami plots."""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import pandas as pd

from ._manhattan_panel import manhattan_spec, render_manhattan_panel
from .backends.base import PlotBackend
from .backends.hover import HoverConfig


@dataclass(frozen=True)
class MiamiRequest:
    """One mirrored Manhattan figure, resolved by the plotter.

    ``top`` and ``bottom`` are prepared against one shared ``GenomeLayout``.
    ``rs_col`` names the id column the annotations index into and is set
    whenever either annotation tuple is non-empty.
    """

    top: pd.DataFrame
    bottom: pd.DataFrame
    hover: Optional[HoverConfig]
    rs_col: Optional[str]
    top_threshold: Optional[float]
    bottom_threshold: Optional[float]
    top_label: Optional[str]
    bottom_label: Optional[str]
    top_annotations: Tuple[str, ...]
    bottom_annotations: Tuple[str, ...]
    highlights: Tuple[Tuple[str, int, int], ...]
    highlight_color: str
    highlight_alpha: float
    figsize: Tuple[float, float]
    title: Optional[str]


def render_miami(backend: PlotBackend, req: MiamiRequest) -> Any:
    """Draw the two mirrored panels, their annotations, and their highlights."""
    fig, axes = backend.create_figure(
        height_ratios=[1.0, 1.0],
        figsize=req.figsize,
        sharex=True,
    )
    top_ax, bottom_ax = axes
    render_manhattan_panel(
        backend,
        top_ax,
        manhattan_spec(
            req.top,
            significance_threshold=req.top_threshold,
            panel_label=req.top_label,
            hover=req.hover,
        ),
    )
    render_manhattan_panel(
        backend,
        bottom_ax,
        manhattan_spec(
            req.bottom,
            significance_threshold=req.bottom_threshold,
            x_label="Chromosome",
            panel_label=req.bottom_label,
            panel_label_y_frac=0.05,
            invert_y=True,
            hover=req.hover,
        ),
    )
    if req.rs_col is not None:
        for ax, frame, ids in (
            (top_ax, req.top, req.top_annotations),
            (bottom_ax, req.bottom, req.bottom_annotations),
        ):
            _add_snp_annotations(backend, ax, frame, req.rs_col, ids)
    if req.highlights:
        offsets = req.top.attrs["layout"].offsets
        for chrom, start, end in req.highlights:
            _add_region_highlight(
                backend,
                [top_ax, bottom_ax],
                str(chrom),
                start,
                end,
                offsets,
                req.highlight_color,
                req.highlight_alpha,
            )
    if req.title:
        backend.set_suptitle(fig, req.title, fontsize=14)
        backend.finalize_layout(fig, top=0.92, hspace=0.05)
    else:
        backend.finalize_layout(fig, hspace=0.05)
    return fig


def _add_snp_annotations(
    backend: PlotBackend,
    ax: Any,
    prepared_df: pd.DataFrame,
    rs_col: str,
    snp_ids: Tuple[str, ...],
) -> None:
    if not snp_ids:
        return
    for _, row in prepared_df[prepared_df[rs_col].isin(snp_ids)].iterrows():
        backend.add_text(
            ax,
            x=row["_cumulative_pos"],
            y=row["neglog10p"],
            text=str(row[rs_col]),
            fontsize=8,
            ha="center",
            va="bottom",
        )


def _add_region_highlight(
    backend: PlotBackend,
    axes: List[Any],
    chrom: str,
    start: int,
    end: int,
    offsets: dict[str, float],
    color: str,
    alpha: float,
) -> None:
    if chrom not in offsets:
        return
    backend.add_region_highlight(
        axes,
        offsets[chrom] + start,
        offsets[chrom] + end,
        color=color,
        alpha=alpha,
    )

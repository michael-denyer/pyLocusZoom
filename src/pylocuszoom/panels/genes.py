"""The gene-track panel: gene bodies, exons, strand arrows, and labels."""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import pandas as pd

from ..backends.base import PlotBackend
from ..colors import GENE_LABEL_COLOR, STRAND_COLORS
from ..config import RegionConfig
from ..gene_track import (
    EXON_HEIGHT,
    GENE_AREA,
    INTRON_HEIGHT,
    ROW_HEIGHT,
    assign_gene_positions,
    compute_arrow_geometry,
    filter_genes_by_region,
)


@dataclass(frozen=True)
class GenePanel:
    """Prepared gene-track panel and its layout.

    ``genes`` and ``exons`` are filtered to the region, ``genes`` is sorted by
    start, and ``rows`` gives the stacking row of each gene in that order.
    """

    genes: pd.DataFrame
    rows: List[int]
    exons: Optional[pd.DataFrame]
    region: RegionConfig
    height: float

    @classmethod
    def from_genes(
        cls,
        genes_df: pd.DataFrame,
        region: RegionConfig,
        exons_df: Optional[pd.DataFrame],
    ) -> "GenePanel":
        """Filter genes to the region and lay them out in non-overlapping rows."""
        genes = filter_genes_by_region(
            genes_df, region.chrom, region.start, region.end
        ).sort_values("start")
        rows = assign_gene_positions(genes, region.start, region.end)
        exons = (
            filter_genes_by_region(exons_df, region.chrom, region.start, region.end)
            if exons_df is not None and not exons_df.empty
            else None
        )
        return cls(
            genes=genes,
            rows=rows,
            exons=exons,
            region=region,
            height=1.0 + max(rows, default=0) * 0.5,
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw gene bodies, exons, strand arrows, and labels for the region."""
        start, end = self.region.start, self.region.end

        backend.set_xlim(ax, start, end)
        backend.set_ylabel(ax, "", fontsize=10)
        backend.hide_yaxis(ax)

        if self.genes.empty:
            backend.set_ylim(ax, 0, 1)
            backend.add_text(
                ax,
                (start + end) / 2,
                0.5,
                "No genes",
                fontsize=9,
                ha="center",
                va="center",
                color="grey",
            )
            return

        max_row = max(self.rows, default=0)
        bottom_margin = EXON_HEIGHT / 2 + 0.02
        backend.set_ylim(ax, -bottom_margin, max_row * ROW_HEIGHT + GENE_AREA + 0.05)

        region_exons = self.exons
        region_width = end - start

        for idx, (_, gene) in enumerate(self.genes.iterrows()):
            gene_start = max(int(gene["start"]), start)
            gene_end = min(int(gene["end"]), end)
            gene_name = gene.get("gene_name", "")

            raw_strand = gene.get("strand")
            strand = raw_strand if raw_strand in ("+", "-") else None
            gene_col = STRAND_COLORS.get(strand, STRAND_COLORS[None])

            y_gene = self.rows[idx] * ROW_HEIGHT + 0.05

            gene_exons = None
            if region_exons is not None and not region_exons.empty and gene_name:
                gene_exons = region_exons[region_exons["gene_name"] == gene_name].copy()

            if gene_exons is not None and not gene_exons.empty:
                _gene_band(
                    backend,
                    ax,
                    (gene_start, gene_end),
                    y_gene,
                    INTRON_HEIGHT,
                    gene_col,
                    1,
                )
                for _, exon in gene_exons.iterrows():
                    span = (max(int(exon["start"]), start), min(int(exon["end"]), end))
                    _gene_band(backend, ax, span, y_gene, EXON_HEIGHT, gene_col, 2)
            else:
                _gene_band(
                    backend,
                    ax,
                    (gene_start, gene_end),
                    y_gene,
                    EXON_HEIGHT,
                    gene_col,
                    2,
                )

            if strand is not None:
                _draw_strand_arrows(
                    backend, ax, strand, gene_start, gene_end, y_gene, region_width
                )

            if gene_name:
                backend.add_text(
                    ax,
                    (gene_start + gene_end) / 2,
                    y_gene + EXON_HEIGHT / 2 + 0.01,
                    gene_name,
                    fontsize=9,
                    ha="center",
                    va="bottom",
                    color=GENE_LABEL_COLOR,
                )


def _gene_band(
    backend: PlotBackend,
    ax: Any,
    span: Tuple[int, int],
    y_gene: float,
    thickness: float,
    color: str,
    zorder: int,
) -> None:
    """Draw one horizontal band of a gene body at ``y_gene``."""
    x0, x1 = span
    backend.add_rectangle(
        ax,
        (x0, y_gene - thickness / 2),
        x1 - x0,
        thickness,
        facecolor=color,
        edgecolor=color,
        linewidth=0.5,
        zorder=zorder,
    )


def _draw_strand_arrows(
    backend: PlotBackend,
    ax: Any,
    strand: str,
    gene_start: int,
    gene_end: int,
    y_gene: float,
    region_width: int,
) -> None:
    """Draw strand direction arrows along a gene body."""
    tips, tri_height, tri_width, arrow_color = compute_arrow_geometry(
        gene_start, gene_end, region_width, strand
    )
    for tip_x in tips:
        base_x = tip_x - tri_width if strand == "+" else tip_x + tri_width
        backend.add_polygon(
            ax,
            [
                [tip_x, y_gene],
                [base_x, y_gene + tri_height],
                [base_x, y_gene - tri_height],
            ],
            facecolor=arrow_color,
            edgecolor=arrow_color,
            linewidth=0.5,
            zorder=5,
        )

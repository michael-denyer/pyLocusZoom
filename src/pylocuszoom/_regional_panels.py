"""Regional panel value types and the drawing for each one.

Every panel resolves its mode, its hover contract, and its layout when it is
built, so the ``draw_*`` functions below take a prepared panel and issue
backend primitives without probing the frame again.  The composer in
``_regional.py`` picks the function; this module owns what it does.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ._plotter_utils import add_significance_line
from .backends.base import (
    PlotBackend,
    SupportsSNPLabels,
)
from .backends.composition import (
    LD_LEGEND_TITLE,
    LegendEntry,
    eqtl_legend_entries,
    finemapping_legend_entries,
    heatmap_highlight_rects,
    ld_legend_entries,
    lower_triangle,
    render_recombination_overlay,
)
from .backends.hover import HoverConfig, HoverDataBuilder
from .colors import (
    LD_HEATMAP_COLORS,
    LEAD_SNP_COLOR,
    LEAD_SNP_HIGHLIGHT_COLOR,
    PIP_LINE_COLOR,
    get_credible_set_color,
    get_eqtl_color,
    get_ld_bin,
    get_ld_color_palette,
)
from .config import ColumnConfig, DisplayConfig, RegionConfig
from .eqtl import prepare_eqtl_for_plotting
from .finemapping import get_credible_sets, prepare_finemapping_for_plotting
from .gene_track import (
    EXON_HEIGHT,
    GENE_AREA,
    INTRON_HEIGHT,
    ROW_HEIGHT,
    STRAND_COLORS,
    assign_gene_positions,
    compute_arrow_geometry,
    filter_genes_by_region,
)
from .logging import logger

REGIONAL_LINE_ALPHA = 0.65
PIP_SCATTER_THRESHOLD = 0.01


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
    height: float
    columns: ColumnConfig
    display: DisplayConfig
    ld_col: Optional[str]
    lead_pos: Optional[int]
    recomb_df: Optional[pd.DataFrame]
    hover: HoverConfig
    panel_label: Optional[str] = None
    add_ld_legend: bool = False


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


@dataclass(frozen=True)
class EqtlPanel:
    """Prepared eQTL panel; ``gene`` is the gene the data was filtered to.

    ``effect_col`` is None when the frame carries no effect sizes, which is
    the single-colour diamond mode.
    """

    data: pd.DataFrame
    height: float
    gene: Optional[str]
    threshold: float
    effect_col: Optional[str]
    hover: HoverConfig

    @classmethod
    def from_frame(
        cls,
        df: pd.DataFrame,
        region: RegionConfig,
        gene: Optional[str],
        threshold: float,
    ) -> "EqtlPanel":
        """Validate, gene- and region-filter, and transform raw eQTL results.

        Raises:
            EQTLValidationError: If required columns are missing, or ``gene``
                is given and the frame has no ``gene`` column.
        """
        data = prepare_eqtl_for_plotting(
            df, gene=gene, chrom=region.chrom, start=region.start, end=region.end
        )
        extra_cols = {
            col: label
            for col, label in (("effect_size", "Effect"), ("gene", "Gene"))
            if col in data.columns
        }
        return cls(
            data=data,
            height=2.0,
            gene=gene,
            threshold=threshold,
            effect_col="effect_size" if "effect_size" in data.columns else None,
            hover=HoverConfig(
                pos_col="pos" if "pos" in data.columns else None,
                p_col="p_value" if "p_value" in data.columns else None,
                extra_cols=extra_cols,
            ),
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
            height=1.0 + max(rows, default=0) * 0.5,
        )


@dataclass(frozen=True)
class HeatmapPanel:
    """Prepared regional LD heatmap panel."""

    matrix: pd.DataFrame
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
            height=height,
            x_positions=x_positions,
            snp_ids=kept_ids,
            metric=metric,
            lead_snp_id=lead_snp_id,
        )


RegionalPanel = Union[
    AssociationPanel,
    FinemappingPanel,
    EqtlPanel,
    GenePanel,
    HeatmapPanel,
]


@dataclass(frozen=True)
class RegionalFigurePlan:
    """Complete ordered plan for one regional figure."""

    chrom: int | str
    start: int
    end: int
    panels: Sequence[RegionalPanel]
    figsize: Tuple[float, float]
    hspace: float = 0.1


def draw_association(
    backend: PlotBackend,
    ax: Any,
    panel: AssociationPanel,
    plan: RegionalFigurePlan,
    genomewide_threshold: float,
) -> None:
    """Draw the association scatter with its axes, overlay, and legends."""
    df = panel.data
    columns = panel.columns
    _draw_association_points(backend, ax, panel)
    add_significance_line(backend, ax, genomewide_threshold, alpha=REGIONAL_LINE_ALPHA)
    backend.set_ylabel(ax, r"$-\log_{10}$ P")
    y_max = df["neglog10p"].max()
    if pd.notna(y_max) and y_max > 0:
        backend.set_ylim(ax, 0, y_max * 1.15)
    backend.set_xlim(ax, plan.start, plan.end)

    if (
        panel.display.snp_labels
        and columns.rs_col in df.columns
        and panel.display.label_top_n > 0
        and not df.empty
        and isinstance(backend, SupportsSNPLabels)
    ):
        backend.add_snp_labels(
            ax,
            df,
            pos_col=columns.pos_col,
            neglog10p_col="neglog10p",
            rs_col=columns.rs_col,
            label_top_n=panel.display.label_top_n,
            adjust=True,
            lead_pos=panel.lead_pos,
            region_span=plan.end - plan.start,
        )

    recomb_df = panel.recomb_df
    if recomb_df is not None and not recomb_df.empty:
        render_recombination_overlay(backend, ax, recomb_df, plan.start, plan.end)

    if panel.panel_label:
        backend.add_panel_label(ax, panel.panel_label)
    if panel.add_ld_legend and panel.ld_col is not None:
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
                colors=palette.get(bin_label, "#BEBEBE"),
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
            colors="#BEBEBE",
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
            "#BEBEBE",
            30,
            0.3,
            2,
        )
    )
    return groups


def draw_finemapping(backend: PlotBackend, ax: Any, panel: FinemappingPanel) -> None:
    """Draw the PIP line, credible-set points, and their legend."""
    data = panel.data
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
        hover_builder = HoverDataBuilder(panel.hover)
        for subset, color, size, linewidth, zorder in _finemapping_groups(panel):
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
        if panel.credible_sets:
            backend.add_legend(
                ax,
                finemapping_legend_entries(panel.credible_sets),
                loc="upper right",
                title="Credible sets",
            )
    backend.set_ylabel(ax, "PIP")
    backend.set_ylim(ax, -0.05, 1.05)


def draw_eqtl(backend: PlotBackend, ax: Any, panel: EqtlPanel) -> None:
    """Draw eQTL points, their legend, and the eQTL significance line."""
    data = panel.data
    if not data.empty:
        hover_builder = HoverDataBuilder(panel.hover)
        if panel.effect_col is not None:
            effect = data[panel.effect_col]
            for subset, marker in (
                (data[effect >= 0], "^"),
                (data[effect < 0], "v"),
            ):
                if not subset.empty:
                    backend.scatter(
                        ax,
                        subset["pos"],
                        subset["neglog10p"],
                        colors=subset[panel.effect_col].apply(get_eqtl_color).tolist(),
                        sizes=50,
                        marker=marker,
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=2,
                        hover_data=hover_builder.build_dataframe(subset),
                    )
            backend.add_legend(
                ax,
                eqtl_legend_entries(),
                loc="upper right",
                title="eQTL effect",
            )
        else:
            label = f"eQTL ({panel.gene})" if panel.gene else "eQTL"
            backend.scatter(
                ax,
                data["pos"],
                data["neglog10p"],
                colors="#FF6B6B",
                sizes=60,
                marker="D",
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
                hover_data=hover_builder.build_dataframe(data),
            )
            backend.add_legend(
                ax,
                [LegendEntry(label, "#FF6B6B", marker="D")],
                loc="upper right",
            )
    backend.set_ylabel(ax, r"$-\log_{10}$ P (eQTL)")
    backend.axhline(
        ax,
        y=-np.log10(panel.threshold),
        color="red",
        linestyle="--",
        linewidth=1,
        alpha=REGIONAL_LINE_ALPHA,
    )


def draw_heatmap(
    backend: PlotBackend, ax: Any, panel: HeatmapPanel, plan: RegionalFigurePlan
) -> None:
    """Draw the lower-triangle LD heatmap and its lead-SNP crosshair."""
    n_snps = len(panel.snp_ids)
    if n_snps < 2:
        logger.debug("Skipping heatmap: fewer than 2 SNPs after filtering")
        return
    mappable = backend.add_heatmap(
        ax,
        data=lower_triangle(panel.matrix.values),
        x_coords=panel.x_positions,
        y_coords=list(range(n_snps)),
        cmap_colors=LD_HEATMAP_COLORS,
        vmin=0.0,
        vmax=1.0,
    )
    backend.add_colorbar(ax, mappable, label="R²" if panel.metric == "r2" else "D'")
    if panel.lead_snp_id is not None and panel.lead_snp_id in panel.snp_ids:
        rects = heatmap_highlight_rects(
            panel.snp_ids.index(panel.lead_snp_id),
            panel.x_positions,
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
    backend.set_xlim(ax, plan.start, plan.end)
    backend.hide_yaxis(ax)


def draw_genes(
    backend: PlotBackend, ax: Any, panel: GenePanel, plan: RegionalFigurePlan
) -> None:
    """Draw gene bodies, exons, strand arrows, and labels for the region."""
    start, end = plan.start, plan.end

    backend.set_xlim(ax, start, end)
    backend.set_ylabel(ax, "", fontsize=10)
    backend.hide_yaxis(ax)

    if panel.genes.empty:
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

    max_row = max(panel.rows, default=0)
    bottom_margin = EXON_HEIGHT / 2 + 0.02
    backend.set_ylim(ax, -bottom_margin, max_row * ROW_HEIGHT + GENE_AREA + 0.05)

    region_exons = panel.exons
    region_width = end - start

    for idx, (_, gene) in enumerate(panel.genes.iterrows()):
        gene_start = max(int(gene["start"]), start)
        gene_end = min(int(gene["end"]), end)
        gene_name = gene.get("gene_name", "")

        raw_strand = gene.get("strand")
        strand = raw_strand if raw_strand in ("+", "-") else None
        gene_col = STRAND_COLORS.get(strand, STRAND_COLORS[None])

        y_gene = panel.rows[idx] * ROW_HEIGHT + 0.05

        gene_exons = None
        if region_exons is not None and not region_exons.empty and gene_name:
            gene_exons = region_exons[region_exons["gene_name"] == gene_name].copy()

        if gene_exons is not None and not gene_exons.empty:
            _gene_band(
                backend, ax, (gene_start, gene_end), y_gene, INTRON_HEIGHT, gene_col, 1
            )
            for _, exon in gene_exons.iterrows():
                span = (max(int(exon["start"]), start), min(int(exon["end"]), end))
                _gene_band(backend, ax, span, y_gene, EXON_HEIGHT, gene_col, 2)
        else:
            _gene_band(
                backend, ax, (gene_start, gene_end), y_gene, EXON_HEIGHT, gene_col, 2
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
                color="#000000",
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

"""Shared regional association composition.

This module owns the policy shared by single and stacked regional plots.  Each
panel type knows how to build itself from raw caller input, and the composer
dispatches on the panel type to draw it, owning the association panel's axes,
labels, optional recombination, and LD legend rules.
"""

from dataclasses import dataclass
from functools import singledispatchmethod
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ._plotter_utils import add_significance_line
from .backends.base import (
    PlotBackend,
    SupportsHeatmap,
    SupportsSecondaryAxis,
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
    get_eqtl_color,
    get_ld_bin,
    get_ld_color_palette,
)
from .config import ColumnConfig, DisplayConfig, RegionConfig
from .eqtl import prepare_eqtl_for_plotting
from .finemapping import plot_finemapping, prepare_finemapping_for_plotting
from .gene_track import (
    assign_gene_positions,
    filter_genes_by_region,
    plot_gene_track_generic,
)
from .logging import logger

REGIONAL_LINE_ALPHA = 0.65


@dataclass(frozen=True)
class AssociationPanel:
    """Prepared association panel and its presentation policy.

    ``data`` already carries ``neglog10p`` and any merged LD column.
    """

    data: pd.DataFrame
    height: float
    columns: ColumnConfig
    display: DisplayConfig
    ld_col: Optional[str]
    lead_pos: Optional[int]
    recomb_df: Optional[pd.DataFrame]
    panel_label: Optional[str] = None
    add_ld_legend: bool = False


@dataclass(frozen=True)
class FinemappingPanel:
    """Prepared fine-mapping panel."""

    data: pd.DataFrame
    height: float
    cs_col: Optional[str]

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
        return cls(data=data, height=1.5, cs_col=cs_col)


@dataclass(frozen=True)
class EqtlPanel:
    """Prepared eQTL panel; ``gene`` is the gene the data was filtered to."""

    data: pd.DataFrame
    height: float
    gene: Optional[str]
    threshold: float

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
        return cls(data=data, height=2.0, gene=gene, threshold=threshold)


@dataclass(frozen=True)
class GenePanel:
    """Prepared gene-track panel; ``data`` is filtered to the region."""

    data: pd.DataFrame
    height: float
    exons_df: Optional[pd.DataFrame]

    @classmethod
    def from_genes(
        cls,
        genes_df: pd.DataFrame,
        region: RegionConfig,
        exons_df: Optional[pd.DataFrame],
    ) -> "GenePanel":
        """Filter genes to the region and size the panel to its stacked rows."""
        genes = filter_genes_by_region(genes_df, region.chrom, region.start, region.end)
        rows = assign_gene_positions(
            genes.sort_values("start"), region.start, region.end
        )
        return cls(
            data=genes, height=1.0 + max(rows, default=0) * 0.5, exons_df=exons_df
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


class RegionalPlotComposer:
    """Compose shared association-panel policy for regional plot modes."""

    def __init__(
        self,
        backend: PlotBackend,
        genomewide_threshold: float,
    ):
        self._backend = backend
        self._genomewide_threshold = genomewide_threshold

    def render(self, plan: RegionalFigurePlan) -> Any:
        """Render every panel in one figure plan and finalize its layout."""
        if not plan.panels:
            raise ValueError("Regional figure plan must contain at least one panel")

        fig, axes = self._backend.create_figure(
            n_panels=len(plan.panels),
            height_ratios=[panel.height for panel in plan.panels],
            figsize=plan.figsize,
            sharex=True,
        )
        for ax, panel in zip(axes, plan.panels):
            self.render_panel(panel, ax, fig, plan)

        self._backend.set_xlabel(axes[-1], f"Chromosome {plan.chrom} (Mb)")
        for ax in axes:
            self._backend.format_xaxis_mb(ax)
        self._backend.finalize_layout(fig, hspace=plan.hspace)
        return fig

    @singledispatchmethod
    def render_panel(
        self, panel: RegionalPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        """Render one panel onto its axes, dispatching on the panel type."""
        raise TypeError(f"No renderer for {type(panel).__name__}")

    @render_panel.register
    def _render_association(
        self, panel: AssociationPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        df = panel.data
        columns = panel.columns
        self.render_association_scatter(
            ax,
            df,
            columns.pos_col,
            panel.ld_col,
            panel.lead_pos,
            columns.rs_col,
            columns.p_col,
        )
        add_significance_line(
            self._backend, ax, self._genomewide_threshold, alpha=REGIONAL_LINE_ALPHA
        )
        self._backend.set_ylabel(ax, r"$-\log_{10}$ P")
        y_max = df["neglog10p"].max()
        if pd.notna(y_max) and y_max > 0:
            self._backend.set_ylim(ax, 0, y_max * 1.15)
        self._backend.set_xlim(ax, plan.start, plan.end)

        if (
            panel.display.snp_labels
            and columns.rs_col in df.columns
            and panel.display.label_top_n > 0
            and not df.empty
            and isinstance(self._backend, SupportsSNPLabels)
        ):
            self._backend.add_snp_labels(
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
        has_recomb = recomb_df is not None and not recomb_df.empty
        if has_recomb and isinstance(self._backend, SupportsSecondaryAxis):
            render_recombination_overlay(
                self._backend, ax, recomb_df, plan.start, plan.end
            )

        if panel.panel_label:
            self._backend.add_panel_label(ax, panel.panel_label)
        if (
            panel.add_ld_legend
            and panel.ld_col is not None
            and panel.ld_col in df.columns
        ):
            self._backend.add_legend(
                ax, ld_legend_entries(), loc="upper right", title=LD_LEGEND_TITLE
            )

    def render_association_scatter(
        self,
        ax: Any,
        df: pd.DataFrame,
        pos_col: str,
        ld_col: Optional[str],
        lead_pos: Optional[int],
        rs_col: str,
        p_col: str,
    ) -> None:
        """Render association points, including LD and lead-SNP styling."""
        hover_config = HoverConfig(
            snp_col=rs_col if rs_col in df.columns else None,
            pos_col=pos_col if pos_col in df.columns else None,
            p_col=p_col if p_col in df.columns else None,
            ld_col=ld_col if ld_col and ld_col in df.columns else None,
        )
        hover_builder = HoverDataBuilder(hover_config)

        if ld_col is not None and ld_col in df.columns:
            df = df.copy()
            df["ld_bin"] = df[ld_col].apply(get_ld_bin)
            df = df.sort_values(ld_col, ascending=True, na_position="first")
            palette = get_ld_color_palette()
            for bin_label in df["ld_bin"].unique():
                bin_data = df[df["ld_bin"] == bin_label]
                self._backend.scatter(
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
            self._backend.scatter(
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

        if lead_pos is not None:
            lead_snp = df[df[pos_col] == lead_pos]
            if not lead_snp.empty:
                self._backend.scatter(
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

    @render_panel.register
    def _render_heatmap(
        self, panel: HeatmapPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        if not isinstance(self._backend, SupportsHeatmap):
            logger.debug(
                "Skipping heatmap: {} does not support heatmaps",
                type(self._backend).__name__,
            )
            return
        n_snps = len(panel.snp_ids)
        if n_snps < 2:
            logger.debug("Skipping heatmap: fewer than 2 SNPs after filtering")
            return
        mappable = self._backend.add_heatmap(
            ax,
            data=lower_triangle(panel.matrix.values),
            x_coords=panel.x_positions,
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
        )
        self._backend.add_colorbar(
            ax, mappable, label="R²" if panel.metric == "r2" else "D'"
        )
        if panel.lead_snp_id is not None and panel.lead_snp_id in panel.snp_ids:
            rects = heatmap_highlight_rects(
                panel.snp_ids.index(panel.lead_snp_id),
                panel.x_positions,
                list(range(n_snps)),
            )
            for x0, y0, width, height in rects:
                self._backend.add_rectangle(
                    ax,
                    (x0, y0),
                    width,
                    height,
                    facecolor=None,
                    edgecolor=LEAD_SNP_HIGHLIGHT_COLOR,
                    linewidth=2,
                    zorder=10,
                )
        self._backend.set_xlim(ax, plan.start, plan.end)
        self._backend.hide_yaxis(ax)

    @render_panel.register
    def _render_genes(
        self, panel: GenePanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        plot_gene_track_generic(
            ax,
            self._backend,
            panel.data,
            plan.chrom,
            plan.start,
            plan.end,
            panel.exons_df,
        )

    @render_panel.register
    def _render_finemapping(
        self, panel: FinemappingPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        fm_data = panel.data
        cs_col = panel.cs_col
        if not fm_data.empty:
            plot_finemapping(
                self._backend,
                ax,
                fm_data,
                pos_col="pos",
                pip_col="pip",
                cs_col=cs_col,
                show_credible_sets=True,
                pip_threshold=0.01,
            )
            credible_sets = []
            if cs_col and cs_col in fm_data.columns:
                credible_sets = sorted(
                    value for value in fm_data[cs_col].dropna().unique() if value != 0
                )
            if credible_sets:
                self._backend.add_legend(
                    ax,
                    finemapping_legend_entries(credible_sets),
                    loc="upper right",
                    title="Credible sets",
                )
        self._backend.set_ylabel(ax, "PIP")
        self._backend.set_ylim(ax, -0.05, 1.05)

    @render_panel.register
    def _render_eqtl(
        self, panel: EqtlPanel, ax: Any, fig: Any, plan: RegionalFigurePlan
    ) -> None:
        eqtl_data = panel.data
        if not eqtl_data.empty:
            extra_cols = {}
            if "effect_size" in eqtl_data.columns:
                extra_cols["effect_size"] = "Effect"
            if "gene" in eqtl_data.columns:
                extra_cols["gene"] = "Gene"
            hover_builder = HoverDataBuilder(
                HoverConfig(
                    pos_col="pos" if "pos" in eqtl_data.columns else None,
                    p_col="p_value" if "p_value" in eqtl_data.columns else None,
                    extra_cols=extra_cols,
                )
            )
            if "effect_size" in eqtl_data.columns:
                pos_effects = eqtl_data[eqtl_data["effect_size"] >= 0]
                neg_effects = eqtl_data[eqtl_data["effect_size"] < 0]
                for subset, marker in ((pos_effects, "^"), (neg_effects, "v")):
                    if not subset.empty:
                        self._backend.scatter(
                            ax,
                            subset["pos"],
                            subset["neglog10p"],
                            colors=subset["effect_size"].apply(get_eqtl_color).tolist(),
                            sizes=50,
                            marker=marker,
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                            hover_data=hover_builder.build_dataframe(subset),
                        )
                self._backend.add_legend(
                    ax,
                    eqtl_legend_entries(),
                    loc="upper right",
                    title="eQTL effect",
                )
            else:
                label = f"eQTL ({panel.gene})" if panel.gene else "eQTL"
                self._backend.scatter(
                    ax,
                    eqtl_data["pos"],
                    eqtl_data["neglog10p"],
                    colors="#FF6B6B",
                    sizes=60,
                    marker="D",
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                    hover_data=hover_builder.build_dataframe(eqtl_data),
                )
                self._backend.add_legend(
                    ax,
                    [LegendEntry(label, "#FF6B6B", marker="D")],
                    loc="upper right",
                )
        self._backend.set_ylabel(ax, r"$-\log_{10}$ P (eQTL)")
        self._backend.axhline(
            ax,
            y=-np.log10(panel.threshold),
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=REGIONAL_LINE_ALPHA,
        )

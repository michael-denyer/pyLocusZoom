"""Shared regional association composition.

This module owns the policy shared by single and stacked regional plots.  The
plotter prepares data and delegates panel composition to this composer, which
owns the association panel's axes, labels, optional recombination, and LD
legend rules.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

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
    ld_legend_entries,
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
from .finemapping import plot_finemapping
from .gene_track import plot_gene_track_generic
from .logging import logger


@dataclass(frozen=True)
class AssociationPanel:
    """Prepared association panel and its presentation policy."""

    data: pd.DataFrame
    height: float
    pos_col: str
    ld_col: Optional[str]
    lead_pos: Optional[int]
    rs_col: Optional[str]
    p_col: Optional[str]
    snp_labels: bool
    label_top_n: int
    genes_df: Optional[pd.DataFrame]
    recomb_df: Optional[pd.DataFrame]
    panel_label: Optional[str] = None
    add_ld_legend: bool = False


@dataclass(frozen=True)
class FinemappingPanel:
    """Prepared fine-mapping panel."""

    data: pd.DataFrame
    height: float
    cs_col: Optional[str]


@dataclass(frozen=True)
class EqtlPanel:
    """Prepared eQTL panel."""

    data: pd.DataFrame
    height: float
    gene_filtered: bool
    gene: Optional[str]
    threshold: float


@dataclass(frozen=True)
class GenePanel:
    """Prepared gene-track panel."""

    data: pd.DataFrame
    height: float
    exons_df: Optional[pd.DataFrame]


@dataclass(frozen=True)
class HeatmapPanel:
    """Prepared regional LD heatmap panel."""

    matrix: pd.DataFrame
    height: float
    x_positions: List[int]
    snp_ids: List[str]
    metric: str
    lead_snp_id: Optional[str]


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
        genomewide_line: float,
    ):
        self._backend = backend
        self._genomewide_line = genomewide_line

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
            if isinstance(panel, AssociationPanel):
                self.render_association_panel(
                    ax,
                    panel.data,
                    pos_col=panel.pos_col,
                    ld_col=panel.ld_col,
                    lead_pos=panel.lead_pos,
                    rs_col=panel.rs_col,
                    p_col=panel.p_col,
                    chrom=plan.chrom,
                    start=plan.start,
                    end=plan.end,
                    snp_labels=panel.snp_labels,
                    label_top_n=panel.label_top_n,
                    genes_df=panel.genes_df,
                    recomb_df=panel.recomb_df,
                    panel_label=panel.panel_label,
                    add_ld_legend=panel.add_ld_legend,
                )
            elif isinstance(panel, FinemappingPanel):
                self.render_finemapping_panel(
                    ax,
                    panel.data,
                    cs_col=panel.cs_col,
                )
            elif isinstance(panel, EqtlPanel):
                self.render_eqtl_panel(
                    ax,
                    panel.data,
                    eqtl_gene_filtered=panel.gene_filtered,
                    eqtl_gene=panel.gene,
                    eqtl_threshold=panel.threshold,
                )
            elif isinstance(panel, GenePanel):
                self.render_gene_panel(
                    ax,
                    panel.data,
                    chrom=plan.chrom,
                    start=plan.start,
                    end=plan.end,
                    exons_df=panel.exons_df,
                )
            elif isinstance(panel, HeatmapPanel):
                self.render_heatmap_panel(
                    ax=ax,
                    fig=fig,
                    ld_matrix=panel.matrix,
                    x_positions=panel.x_positions,
                    snp_ids=panel.snp_ids,
                    metric=panel.metric,
                    lead_snp_id=panel.lead_snp_id,
                    start=plan.start,
                    end=plan.end,
                )

        self._backend.set_xlabel(axes[-1], f"Chromosome {plan.chrom} (Mb)")
        for ax in axes:
            self._backend.format_xaxis_mb(ax)
        self._backend.finalize_layout(fig, hspace=plan.hspace)
        return fig

    def render_association_panel(
        self,
        ax: Any,
        df: pd.DataFrame,
        *,
        pos_col: str,
        ld_col: Optional[str],
        lead_pos: Optional[int],
        rs_col: Optional[str],
        p_col: Optional[str],
        chrom: int | str,
        start: int,
        end: int,
        snp_labels: bool,
        label_top_n: int,
        genes_df: Optional[pd.DataFrame],
        recomb_df: Optional[pd.DataFrame],
        panel_label: Optional[str] = None,
        add_ld_legend: bool = False,
    ) -> None:
        """Render one association panel with shared regional policy."""
        self.render_association_scatter(
            ax, df, pos_col, ld_col, lead_pos, rs_col, p_col
        )
        self._backend.axhline(
            ax,
            y=self._genomewide_line,
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=0.65,
            zorder=1,
        )
        self._backend.set_ylabel(ax, r"$-\log_{10}$ P")
        y_max = df["neglog10p"].max()
        if pd.notna(y_max) and y_max > 0:
            self._backend.set_ylim(ax, 0, y_max * 1.15)
        self._backend.set_xlim(ax, start, end)

        if (
            snp_labels
            and rs_col is not None
            and rs_col in df.columns
            and label_top_n > 0
            and not df.empty
            and isinstance(self._backend, SupportsSNPLabels)
        ):
            self._backend.add_snp_labels(
                ax,
                df,
                pos_col=pos_col,
                neglog10p_col="neglog10p",
                rs_col=rs_col,
                label_top_n=label_top_n,
                genes_df=genes_df,
                chrom=chrom,
                adjust=True,
                lead_pos=lead_pos,
                region_span=end - start,
            )

        has_recomb = recomb_df is not None and not recomb_df.empty
        if has_recomb and isinstance(self._backend, SupportsSecondaryAxis):
            render_recombination_overlay(self._backend, ax, recomb_df, start, end)

        if panel_label:
            self._backend.add_panel_label(ax, panel_label)
        if add_ld_legend and ld_col is not None and ld_col in df.columns:
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
        rs_col: Optional[str] = None,
        p_col: Optional[str] = None,
    ) -> None:
        """Render association points, including LD and lead-SNP styling."""
        hover_config = HoverConfig(
            snp_col=rs_col if rs_col and rs_col in df.columns else None,
            pos_col=pos_col if pos_col in df.columns else None,
            p_col=p_col if p_col and p_col in df.columns else None,
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

    def render_heatmap_panel(
        self,
        *,
        ax: Any,
        fig: Any,
        ld_matrix: pd.DataFrame,
        x_positions: List[int],
        snp_ids: List[str],
        metric: str,
        lead_snp_id: Optional[str],
        start: int,
        end: int,
    ) -> None:
        """Render an LD heatmap panel in genomic coordinates."""
        if not isinstance(self._backend, SupportsHeatmap):
            logger.debug(
                "Skipping heatmap: {} does not support heatmaps",
                type(self._backend).__name__,
            )
            return
        n_snps = len(snp_ids)
        if n_snps < 2:
            logger.debug("Skipping heatmap: fewer than 2 SNPs after filtering")
            return
        mappable = self._backend.add_heatmap(
            ax,
            data=ld_matrix.values,
            x_coords=x_positions,
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
            mask_upper=True,
        )
        self._backend.add_colorbar(ax, mappable, label="R²" if metric == "r2" else "D'")
        if lead_snp_id is not None and lead_snp_id in snp_ids:
            self._backend.highlight_heatmap_snp(
                ax,
                fig,
                snp_ids.index(lead_snp_id),
                n_snps,
                color=LEAD_SNP_HIGHLIGHT_COLOR,
                linewidth=2,
            )
        self._backend.set_xlim(ax, start, end)
        self._backend.hide_yaxis(ax)

    def render_gene_panel(
        self,
        ax: Any,
        genes_df: pd.DataFrame,
        *,
        chrom: int | str,
        start: int,
        end: int,
        exons_df: Optional[pd.DataFrame],
    ) -> None:
        """Render a gene panel and apply its shared axis policy."""
        plot_gene_track_generic(
            ax, self._backend, genes_df, chrom, start, end, exons_df
        )

    def render_finemapping_panel(
        self,
        ax: Any,
        fm_data: pd.DataFrame,
        *,
        cs_col: Optional[str],
    ) -> None:
        """Render a fine-mapping panel and shared panel policy."""
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

    def render_eqtl_panel(
        self,
        ax: Any,
        eqtl_data: pd.DataFrame,
        *,
        eqtl_gene_filtered: bool,
        eqtl_gene: Optional[str],
        eqtl_threshold: float,
    ) -> None:
        """Render prepared eQTL points and their semantic legend policy."""
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
                label = f"eQTL ({eqtl_gene})" if eqtl_gene_filtered else "eQTL"
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
            y=-np.log10(eqtl_threshold),
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=0.65,
        )

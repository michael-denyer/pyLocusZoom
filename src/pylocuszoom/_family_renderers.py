"""Semantic renderers for plotter families beyond Manhattan/QQ.

These classes accept prepared data and figure intent.  Plotter classes keep
validation and domain preparation; backend primitives stay behind these
rendering seams.
"""

from typing import Any, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from ._plotter_utils import (
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
    add_significance_line,
)
from .backends.base import PlotBackend
from .backends.hover import HoverConfig, HoverDataBuilder
from .colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    LD_BINS,
    LD_HEATMAP_COLORS,
    LD_NA_COLOR,
    LEAD_SNP_COLOR,
    LEAD_SNP_HIGHLIGHT_COLOR,
    SECONDARY_HIGHLIGHT_COLOR,
    get_phewas_category_palette,
)


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
        self._backend.add_region_highlight(
            fig,
            [top_ax, bottom_ax],
            offsets[chrom] + start,
            offsets[chrom] + end,
            color=color,
            alpha=alpha,
        )


class StatsRenderer:
    """Render prepared PheWAS and forest figure intents."""

    def __init__(self, backend: PlotBackend):
        self._backend = backend

    def render_phewas(
        self,
        df: pd.DataFrame,
        *,
        variant_id: str,
        phenotype_col: str,
        p_col: str,
        category_col: str,
        effect_col: Optional[str],
        significance_threshold: float,
        figsize: Tuple[float, float],
    ) -> Any:
        if category_col in df.columns:
            df = df.sort_values([category_col, p_col])
            categories = df[category_col].unique().tolist()
            palette = get_phewas_category_palette(categories)
        else:
            df = df.sort_values(p_col)
            categories, palette = [], {}
        df = df.copy()
        df["y_pos"] = range(len(df))
        fig, axes = self._backend.create_figure(1, [1.0], figsize=figsize)
        ax = axes[0]
        groups = categories or [None]
        for cat in groups:
            if cat is None:
                cat_data, color = df, "#4169E1"
            elif pd.isna(cat):
                cat_data, color = df[df[category_col].isna()], palette[cat]
            else:
                cat_data, color = df[df[category_col] == cat], palette[cat]
            if cat_data.empty:
                continue
            if effect_col and effect_col in cat_data.columns:
                effect_vals = cat_data[effect_col]
                subsets = [
                    (effect_vals.isna(), "o"),
                    (effect_vals >= 0, "^"),
                    (effect_vals < 0, "v"),
                ]
                for mask, marker in subsets:
                    subset = cat_data[mask]
                    if not subset.empty:
                        self._backend.scatter(
                            ax,
                            subset["neglog10p"],
                            subset["y_pos"],
                            colors=color,
                            sizes=60,
                            marker=marker,
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
            else:
                self._backend.scatter(
                    ax,
                    cat_data["neglog10p"],
                    cat_data["y_pos"],
                    colors=color,
                    sizes=60,
                    marker="o",
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                )
        self._backend.axvline(
            ax,
            x=-np.log10(significance_threshold),
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )
        self._backend.set_xlabel(ax, r"$-\log_{10}$ P")
        self._backend.set_ylabel(ax, "Phenotype")
        self._backend.set_ylim(ax, -0.5, len(df) - 0.5)
        self._backend.set_yticks(
            ax, df["y_pos"].tolist(), df[phenotype_col].tolist(), fontsize=8
        )
        self._backend.set_title(ax, f"PheWAS: {variant_id}")
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)
        return fig

    def render_forest(
        self,
        df: pd.DataFrame,
        *,
        variant_id: str,
        study_col: str,
        effect_col: str,
        ci_lower_col: str,
        ci_upper_col: str,
        weight_col: Optional[str],
        null_value: float,
        effect_label: str,
        figsize: Tuple[float, float],
    ) -> Any:
        df = df.copy()
        df["y_pos"] = range(len(df) - 1, -1, -1)
        if weight_col and weight_col in df.columns:
            weights = df[weight_col]
            weight_range = weights.max() - weights.min()
            sizes = (
                40 + (weights - weights.min()) / weight_range * 160
                if weight_range > 0
                else 120
            )
        else:
            sizes = 80
        fig, axes = self._backend.create_figure(1, [1.0], figsize=figsize)
        ax = axes[0]
        self._backend.errorbar_h(
            ax,
            x=df[effect_col],
            y=df["y_pos"],
            xerr_lower=df[effect_col] - df[ci_lower_col],
            xerr_upper=df[ci_upper_col] - df[effect_col],
            color="black",
            linewidth=1.5,
            capsize=3,
            zorder=2,
        )
        self._backend.scatter(
            ax,
            df[effect_col],
            df["y_pos"],
            colors="#4169E1",
            sizes=sizes,
            marker="s",
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )
        self._backend.axvline(
            ax, x=null_value, color="grey", linestyle="--", linewidth=1, alpha=0.7
        )
        self._backend.set_xlabel(ax, effect_label)
        self._backend.set_ylim(ax, -0.5, len(df) - 0.5)
        x_min = min(df[ci_lower_col].min(), null_value)
        x_max = max(df[ci_upper_col].max(), null_value)
        padding = (x_max - x_min) * 0.1
        self._backend.set_xlim(ax, x_min - padding, x_max + padding)
        self._backend.set_yticks(
            ax, df["y_pos"].tolist(), df[study_col].tolist(), fontsize=10
        )
        self._backend.set_title(ax, f"Forest Plot: {variant_id}")
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)
        return fig


class LDHeatmapRenderer:
    """Render a prepared LD matrix intent."""

    def __init__(self, backend: PlotBackend):
        self._backend = backend

    def render(
        self,
        data: np.ndarray,
        snp_ids: List[str],
        *,
        lead_idx: Optional[int],
        highlight_indices: List[int],
        metric: str,
        figsize: Tuple[float, float],
        title: Optional[str],
        show_colorbar: bool,
    ) -> Any:
        n_snps = len(snp_ids)
        fig, axes = self._backend.create_figure(1, [1.0], figsize=figsize, sharex=False)
        ax = axes[0]
        mappable = self._backend.add_heatmap(
            ax,
            data=data,
            x_coords=list(range(n_snps)),
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
            mask_upper=True,
        )
        if show_colorbar:
            self._backend.add_colorbar(
                ax, mappable, label="R²" if metric == "r2" else "D'"
            )
        if lead_idx is not None:
            self._highlight(ax, fig, lead_idx, n_snps, LEAD_SNP_HIGHLIGHT_COLOR)
        for idx in highlight_indices:
            self._highlight(ax, fig, idx, n_snps, SECONDARY_HIGHLIGHT_COLOR)
        ticks = list(range(n_snps))
        self._backend.set_xticks(ax, ticks, snp_ids, rotation=90)
        self._backend.set_yticks(ax, ticks, snp_ids)
        if title:
            self._backend.set_title(ax, title)
        self._backend.finalize_layout(fig)
        return fig

    def _highlight(self, ax: Any, fig: Any, idx: int, n_snps: int, color: str) -> None:
        self._backend.highlight_heatmap_snp(
            ax, fig, idx, n_snps, color=color, linewidth=2
        )


class ColocRenderer:
    """Render a prepared colocalization scatter intent."""

    def __init__(self, backend: PlotBackend):
        self._backend = backend

    def render(
        self,
        merged: pd.DataFrame,
        *,
        merged_rs_col: Optional[str],
        ld_col_merged: Optional[str],
        lead_idx: Optional[Any],
        gwas_threshold: float,
        eqtl_threshold: float,
        show_correlation: bool,
        color_by_effect: bool,
        h4_posterior: Optional[float],
        title: Optional[str],
        figsize: Tuple[float, float],
    ) -> Any:
        fig, axes = self._backend.create_figure(1, [1.0], figsize=figsize)
        ax = axes[0]
        if lead_idx is not None:
            lead_row, other_rows = merged.loc[[lead_idx]], merged.drop(lead_idx)
        else:
            lead_row, other_rows = pd.DataFrame(), merged
        if not other_rows.empty:
            self._backend.scatter(
                ax,
                other_rows["neglog10_gwas"],
                other_rows["neglog10_eqtl"],
                colors=other_rows["color"].tolist(),
                sizes=60,
                marker="o",
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
            )
        if not lead_row.empty:
            self._backend.scatter(
                ax,
                lead_row["neglog10_gwas"],
                lead_row["neglog10_eqtl"],
                colors=LEAD_SNP_COLOR,
                sizes=100,
                marker="D",
                edgecolor="black",
                linewidth=0.5,
                zorder=5,
            )
            if merged_rs_col is not None:
                self._backend.add_text(
                    ax,
                    lead_row["neglog10_gwas"].values[0],
                    lead_row["neglog10_eqtl"].values[0] + 0.5,
                    str(lead_row[merged_rs_col].values[0]),
                    fontsize=9,
                    ha="center",
                    va="bottom",
                )
        self._backend.axvline(
            ax,
            x=-np.log10(gwas_threshold),
            color="grey",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )
        self._backend.axhline(
            ax,
            y=-np.log10(eqtl_threshold),
            color="grey",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )
        x_min, x_max = merged["neglog10_gwas"].min(), merged["neglog10_gwas"].max()
        y_min, y_max = merged["neglog10_eqtl"].min(), merged["neglog10_eqtl"].max()
        x_range, y_range = x_max - x_min, y_max - y_min
        if show_correlation and len(merged) >= 3:
            r, p = stats.pearsonr(merged["neglog10_gwas"], merged["neglog10_eqtl"])
            self._backend.add_text(
                ax,
                x_min + 0.05 * x_range,
                y_max - 0.05 * y_range,
                f"r = {r:.3f}\n{'p < 0.001' if p < 0.001 else f'p = {p:.3f}'}",
                fontsize=10,
                ha="left",
                va="top",
            )
        if h4_posterior is not None:
            self._backend.add_text(
                ax,
                x_max - 0.05 * x_range,
                y_min + 0.05 * y_range,
                f"H4 PP = {h4_posterior:.3f}",
                fontsize=10,
                ha="right",
                va="bottom",
            )
        self._backend.set_xlabel(ax, r"GWAS $-\log_{10}$ P")
        self._backend.set_ylabel(ax, r"eQTL $-\log_{10}$ P")
        if title:
            self._backend.set_title(ax, title)
        self._backend.hide_spines(ax, ["top", "right"])
        if color_by_effect:
            self._backend.add_effect_legend(
                ax,
                [
                    (0.0, "Same direction", EFFECT_CONGRUENT_COLOR),
                    (0.0, "Opposite direction", EFFECT_INCONGRUENT_COLOR),
                    (0.0, "Missing effect", LD_NA_COLOR),
                ],
            )
        elif ld_col_merged is not None:
            self._backend.add_ld_legend(ax, LD_BINS, LEAD_SNP_COLOR)
        self._backend.finalize_layout(fig)
        return fig

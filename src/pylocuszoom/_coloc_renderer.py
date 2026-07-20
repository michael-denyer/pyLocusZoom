"""Semantic renderer for colocalization plots."""

from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from .backends.base import PlotBackend
from .colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    LD_BINS,
    LD_NA_COLOR,
    LEAD_SNP_COLOR,
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

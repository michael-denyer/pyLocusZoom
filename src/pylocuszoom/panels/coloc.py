"""The colocalization panel, which draws itself."""

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd
from scipy import stats

from .._plotter_utils import add_significance_line
from ..backends.base import PlotBackend
from ..backends.composition import (
    LD_LEGEND_TITLE,
    effect_legend_entries,
    ld_legend_entries,
)
from ..colors import LEAD_SNP_COLOR
from ..config import ColocConfig


@dataclass(frozen=True)
class ColocPanel:
    """One colocalization scatter, resolved by the plotter.

    ``rs_col`` and ``ld_col`` are the merged frame's names for the columns
    ``config`` names in the input frames, or None when they are absent.
    """

    merged: pd.DataFrame
    config: ColocConfig
    rs_col: Optional[str]
    ld_col: Optional[str]
    lead_idx: Optional[Any]
    title: Optional[str]

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the GWAS-versus-eQTL scatter, its thresholds, and its legend."""
        merged = self.merged
        config = self.config
        if self.lead_idx is not None:
            lead_row, other_rows = (
                merged.loc[[self.lead_idx]],
                merged.drop(self.lead_idx),
            )
        else:
            lead_row, other_rows = pd.DataFrame(), merged
        if not other_rows.empty:
            backend.scatter(
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
            backend.scatter(
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
            if self.rs_col is not None:
                backend.add_text(
                    ax,
                    lead_row["neglog10_gwas"].values[0],
                    lead_row["neglog10_eqtl"].values[0] + 0.5,
                    str(lead_row[self.rs_col].values[0]),
                    fontsize=9,
                    ha="center",
                    va="bottom",
                )
        add_significance_line(
            backend, ax, config.gwas_threshold, axis="x", color="grey", alpha=0.7
        )
        add_significance_line(
            backend, ax, config.eqtl_threshold, axis="y", color="grey", alpha=0.7
        )
        x_min, x_max = merged["neglog10_gwas"].min(), merged["neglog10_gwas"].max()
        y_min, y_max = merged["neglog10_eqtl"].min(), merged["neglog10_eqtl"].max()
        x_range, y_range = x_max - x_min, y_max - y_min
        if config.show_correlation and len(merged) >= 3:
            r, p = stats.pearsonr(merged["neglog10_gwas"], merged["neglog10_eqtl"])
            backend.add_text(
                ax,
                x_min + 0.05 * x_range,
                y_max - 0.05 * y_range,
                f"r = {r:.3f}\n{'p < 0.001' if p < 0.001 else f'p = {p:.3f}'}",
                fontsize=10,
                ha="left",
                va="top",
            )
        if config.h4_posterior is not None:
            backend.add_text(
                ax,
                x_max - 0.05 * x_range,
                y_min + 0.05 * y_range,
                f"H4 PP = {config.h4_posterior:.3f}",
                fontsize=10,
                ha="right",
                va="bottom",
            )
        backend.set_xlabel(ax, r"GWAS $-\log_{10}$ P")
        backend.set_ylabel(ax, r"eQTL $-\log_{10}$ P")
        if self.title:
            backend.set_title(ax, self.title)
        if config.color_by_effect:
            backend.add_legend(
                ax, effect_legend_entries(), loc="upper right", title="Effect"
            )
        elif self.ld_col is not None:
            backend.add_legend(
                ax, ld_legend_entries(), loc="upper right", title=LD_LEGEND_TITLE
            )

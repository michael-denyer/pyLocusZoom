"""The colocalization panel: the merged GWAS-eQTL frame, which draws itself."""

from dataclasses import dataclass
from typing import Any, Literal, Optional, Tuple

import pandas as pd
from scipy import stats

from .._data import prepare_pvalue_data
from .._plotter_utils import add_significance_line
from ..backends.base import PlotBackend
from ..backends.composition import (
    LD_LEGEND_TITLE,
    effect_legend_entries,
    ld_legend_entries,
)
from ..colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    LD_NA_COLOR,
    LEAD_SNP_COLOR,
    get_ld_color,
)
from ..config import ColocConfig
from ..exceptions import ValidationError
from ..schemas import Canonical, coloc_plot_spec
from ..validation import check, resolve_column


def _get_effect_agreement_color(gwas_effect: float, eqtl_effect: float) -> str:
    """Get color based on effect direction agreement.

    Args:
        gwas_effect: GWAS effect size (beta coefficient).
        eqtl_effect: eQTL effect size (beta coefficient).

    Returns:
        Hex color code: green for same direction, red for opposite.
    """
    if pd.isna(gwas_effect) or pd.isna(eqtl_effect):
        return LD_NA_COLOR
    same_direction = (gwas_effect > 0) == (eqtl_effect > 0)
    return EFFECT_CONGRUENT_COLOR if same_direction else EFFECT_INCONGRUENT_COLOR


def _project_coloc_input(
    df: pd.DataFrame,
    *,
    name: str,
    pos_col: str,
    p_col: str,
    effect_col: Optional[str] = None,
    rs_col: Optional[str] = None,
    ld_col: Optional[str] = None,
) -> pd.DataFrame:
    """Select roles from their declared source before any merge can rename them."""
    roles = {"pos": pos_col, f"p_{name}": p_col}
    for role, source, field, default in (
        (f"{name}_effect", effect_col, f"{name}_effect_col", None),
        ("ld", ld_col, "ld_col", None),
        ("rs", rs_col, "rs_col", Canonical.RS),
    ):
        resolved = resolve_column(
            df,
            source,
            parameter=field,
            optional_default=default,
            frame=f"the {name.upper()} data",
        )
        if resolved is not None:
            roles[role] = resolved
    projected = pd.DataFrame({role: df[source] for role, source in roles.items()})
    return prepare_pvalue_data(
        projected, f"p_{name}", "coloc", out_col=f"neglog10_{name}"
    )


def _merge_and_transform(
    gwas_df: pd.DataFrame,
    eqtl_df: pd.DataFrame,
    config: ColocConfig,
) -> pd.DataFrame:
    """Project source-owned roles, then merge and colour the accepted rows."""
    gwas_effect, eqtl_effect = (
        (config.gwas_effect_col, config.eqtl_effect_col)
        if config.color_by_effect
        else (None, None)
    )
    gwas = _project_coloc_input(
        gwas_df,
        name="gwas",
        pos_col=config.pos_col,
        p_col=config.gwas_p_col,
        effect_col=gwas_effect,
        rs_col=config.rs_col,
        ld_col=config.ld_col,
    )
    eqtl = _project_coloc_input(
        eqtl_df,
        name="eqtl",
        pos_col=config.pos_col,
        p_col=config.eqtl_p_col,
        effect_col=eqtl_effect,
    )
    merged = pd.merge(gwas, eqtl, on="pos", how="inner")
    if merged.empty:
        raise ValidationError(
            "No overlapping positions between GWAS and eQTL DataFrames"
        )
    if config.color_by_effect:
        merged["color"] = merged.apply(
            lambda row: _get_effect_agreement_color(
                row["gwas_effect"], row["eqtl_effect"]
            ),
            axis=1,
        )
    elif "ld" in merged:
        merged["color"] = merged["ld"].apply(get_ld_color)
    else:
        merged["color"] = LD_NA_COLOR
    return merged


def _resolve_lead_idx(merged: pd.DataFrame, config: ColocConfig) -> Optional[Any]:
    """Find the row to draw as the lead variant.

    A named ``lead_snp`` wins. Otherwise a lead is auto-selected by highest
    combined signal, but only when LD colouring is in play: without it there
    is no gradient for the lead to anchor.

    Args:
        merged: Output of ``_merge_and_transform``.
        config: Validated plot configuration.

    Returns:
        Index label of the lead row, or None to draw no lead marker.

    Raises:
        ValidationError: If ``lead_snp`` is named but the merged frame has no SNP
            ID column or no row matching it.
    """
    if config.lead_snp is not None:
        if "rs" not in merged:
            raise ValidationError(
                f"lead_snp '{config.lead_snp}' specified but rs_col not found"
            )
        matches = merged[merged["rs"] == config.lead_snp]
        if len(matches) == 0:
            raise ValidationError(
                f"lead_snp '{config.lead_snp}' not found in merged data"
            )
        return matches.index[0]
    if "ld" in merged:
        combined = merged["neglog10_gwas"] + merged["neglog10_eqtl"]
        return combined.idxmax()
    return None


@dataclass(frozen=True)
class ColocPanel:
    """One colocalization scatter, resolved when it is built.

    ``merged`` carries ``neglog10_gwas``, ``neglog10_eqtl`` and ``color``.
    ``lead_idx`` labels its lead row, and ``lead_label`` is the SNP id written
    beside it, or None when the frame has no ids. ``correlation`` is Pearson's
    ``(r, p)`` over the two p-value columns, or None when it is not shown.
    """

    merged: pd.DataFrame
    lead_idx: Optional[Any]
    lead_label: Optional[str]
    legend: Literal["effect", "ld", None]
    gwas_threshold: Optional[float]
    eqtl_threshold: Optional[float]
    correlation: Optional[Tuple[float, float]]
    h4_posterior: Optional[float]
    title: Optional[str]

    @classmethod
    def from_frames(
        cls,
        gwas_df: pd.DataFrame,
        eqtl_df: pd.DataFrame,
        config: ColocConfig,
        *,
        gwas_threshold: Optional[float],
        eqtl_threshold: Optional[float],
        title: Optional[str],
    ) -> "ColocPanel":
        """Validate both frames, merge them on position, and resolve the lead.

        Raises:
            ValidationError: If a required or named column is missing or
                invalid, no position is in both frames, or ``lead_snp`` is not
                found.
        """
        for frame, name, p_col in (
            (gwas_df, "GWAS DataFrame", config.gwas_p_col),
            (eqtl_df, "eQTL DataFrame", config.eqtl_p_col),
        ):
            check(frame, coloc_plot_spec(name, config.pos_col, p_col))
        merged = _merge_and_transform(gwas_df, eqtl_df, config)
        lead_idx = _resolve_lead_idx(merged, config)
        correlation = None
        if config.show_correlation and len(merged) >= 3:
            r, p = stats.pearsonr(merged["neglog10_gwas"], merged["neglog10_eqtl"])
            correlation = (float(r), float(p))
        if config.color_by_effect:
            legend = "effect"
        else:
            legend = "ld" if "ld" in merged else None
        has_label = lead_idx is not None and "rs" in merged
        return cls(
            merged=merged,
            lead_idx=lead_idx,
            lead_label=str(merged.at[lead_idx, "rs"]) if has_label else None,
            legend=legend,
            gwas_threshold=gwas_threshold,
            eqtl_threshold=eqtl_threshold,
            correlation=correlation,
            h4_posterior=config.h4_posterior,
            title=title,
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the GWAS-versus-eQTL scatter, its thresholds, and its legend."""
        merged = self.merged
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
            if self.lead_label is not None:
                backend.add_text(
                    ax,
                    lead_row["neglog10_gwas"].values[0],
                    lead_row["neglog10_eqtl"].values[0] + 0.5,
                    self.lead_label,
                    fontsize=9,
                    ha="center",
                    va="bottom",
                )
        add_significance_line(
            backend, ax, self.gwas_threshold, axis="x", color="grey", alpha=0.7
        )
        add_significance_line(
            backend, ax, self.eqtl_threshold, axis="y", color="grey", alpha=0.7
        )
        x_min, x_max = merged["neglog10_gwas"].min(), merged["neglog10_gwas"].max()
        y_min, y_max = merged["neglog10_eqtl"].min(), merged["neglog10_eqtl"].max()
        x_range, y_range = x_max - x_min, y_max - y_min
        if self.correlation is not None:
            r, p = self.correlation
            backend.add_text(
                ax,
                x_min + 0.05 * x_range,
                y_max - 0.05 * y_range,
                f"r = {r:.3f}\n{'p < 0.001' if p < 0.001 else f'p = {p:.3f}'}",
                fontsize=10,
                ha="left",
                va="top",
            )
        if self.h4_posterior is not None:
            backend.add_text(
                ax,
                x_max - 0.05 * x_range,
                y_min + 0.05 * y_range,
                f"H4 PP = {self.h4_posterior:.3f}",
                fontsize=10,
                ha="right",
                va="bottom",
            )
        backend.set_xlabel(ax, r"GWAS $-\log_{10}$ P")
        backend.set_ylabel(ax, r"eQTL $-\log_{10}$ P")
        if self.title:
            backend.set_title(ax, self.title)
        if self.legend == "effect":
            backend.add_legend(ax, effect_legend_entries(), title="Effect")
        elif self.legend == "ld":
            backend.add_legend(ax, ld_legend_entries(), title=LD_LEGEND_TITLE)

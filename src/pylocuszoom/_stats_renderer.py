"""Semantic renderers for statistical plot families."""

from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd

from .backends.base import PlotBackend, SupportsBarCharts
from .colors import get_phewas_category_palette


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
        significance_threshold: Optional[float],
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
        fig, axes = self._backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=figsize
        )
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
        if significance_threshold is not None:
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
        """Draw a forest plot of per-study effects and confidence intervals.

        Raises:
            TypeError: If the backend does not implement ``SupportsBarCharts``.
                The bars are the whole figure here, so there is nothing to
                degrade to.
        """
        if not isinstance(self._backend, SupportsBarCharts):
            raise TypeError(
                f"{type(self._backend).__name__} does not support bar charts. "
                "A forest plot needs a backend implementing SupportsBarCharts "
                "(hbar, errorbar_h)."
            )
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
        fig, axes = self._backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=figsize
        )
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
        self._backend.finalize_layout(fig)
        return fig

"""PheWAS and forest panels, each of which draws itself."""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple, Union

import pandas as pd

from ._plotter_utils import add_significance_line
from .backends.base import PlotBackend
from .colors import get_phewas_category_palette

UNCATEGORISED = "Uncategorised"
UNCATEGORISED_COLOR = "#4169E1"


def _phewas_groups(
    df: pd.DataFrame, category_col: str, p_col: str
) -> Tuple[pd.DataFrame, dict[str, str]]:
    """Assign every PheWAS row a drawing group and give each group a colour.

    Args:
        df: Validated PheWAS results.
        category_col: Column naming each phenotype's category, absent from the
            frame when the caller supplied none.
        p_col: Column name for p-value.

    Returns:
        A sorted copy carrying ``_group``, and one hex colour per group.
    """
    if category_col in df.columns:
        df = df.sort_values([category_col, p_col]).copy()
        df["_group"] = df[category_col].fillna(UNCATEGORISED)
        return df, get_phewas_category_palette(df["_group"].unique().tolist())
    df = df.sort_values(p_col).copy()
    df["_group"] = UNCATEGORISED
    return df, {UNCATEGORISED: UNCATEGORISED_COLOR}


def _effect_subsets(
    data: pd.DataFrame, effect_col: Optional[str]
) -> List[Tuple[pd.Series, str]]:
    """Split one PheWAS group into the marker shapes it is drawn with.

    Args:
        data: Rows of one category.
        effect_col: Column holding effect direction, or None to draw one shape.

    Returns:
        One (row mask, marker) pair per shape.
    """
    if effect_col and effect_col in data.columns:
        effects = data[effect_col]
        return [(effects.isna(), "o"), (effects >= 0, "^"), (effects < 0, "v")]
    return [(pd.Series(True, index=data.index), "o")]


@dataclass(frozen=True)
class PhewasPanel:
    """Prepared PheWAS panel.

    ``data`` is in drawing order and carries ``_group`` and ``y_pos``;
    ``palette`` gives one colour per group.
    """

    data: pd.DataFrame
    palette: dict[str, str]
    phenotype_col: str
    effect_col: Optional[str]
    variant_id: str
    significance_threshold: Optional[float]

    @classmethod
    def from_frame(
        cls,
        df: pd.DataFrame,
        *,
        variant_id: str,
        phenotype_col: str,
        p_col: str,
        category_col: str,
        effect_col: Optional[str],
        significance_threshold: Optional[float],
    ) -> "PhewasPanel":
        """Group a validated, p-value-prepared PheWAS frame for drawing."""
        data, palette = _phewas_groups(df, category_col, p_col)
        data["y_pos"] = range(len(data))
        return cls(
            data=data,
            palette=palette,
            phenotype_col=phenotype_col,
            effect_col=effect_col,
            variant_id=variant_id,
            significance_threshold=significance_threshold,
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw one scatter per category, the significance line, and the axes."""
        df = self.data
        for group, group_data in df.groupby("_group", sort=False):
            for mask, marker in _effect_subsets(group_data, self.effect_col):
                subset = group_data[mask]
                if subset.empty:
                    continue
                backend.scatter(
                    ax,
                    subset["neglog10p"],
                    subset["y_pos"],
                    colors=self.palette[group],
                    sizes=60,
                    marker=marker,
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                )
        add_significance_line(
            backend, ax, self.significance_threshold, axis="x", alpha=0.7
        )
        backend.set_xlabel(ax, r"$-\log_{10}$ P")
        backend.set_ylabel(ax, "Phenotype")
        backend.set_ylim(ax, -0.5, len(df) - 0.5)
        backend.set_yticks(
            ax, df["y_pos"].tolist(), df[self.phenotype_col].tolist(), fontsize=8
        )
        backend.set_title(ax, f"PheWAS: {self.variant_id}")


@dataclass(frozen=True)
class ForestPanel:
    """Prepared forest panel.

    ``data`` carries ``y_pos`` with the first study on top; ``sizes`` is one
    marker size per row, or one size for every row when there are no weights.
    """

    data: pd.DataFrame
    sizes: Union[float, pd.Series]
    study_col: str
    effect_col: str
    ci_lower_col: str
    ci_upper_col: str
    null_value: float
    effect_label: str
    variant_id: str

    @classmethod
    def from_frame(
        cls,
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
    ) -> "ForestPanel":
        """Lay out a validated forest frame and size its markers by weight."""
        data = df.copy()
        data["y_pos"] = range(len(data) - 1, -1, -1)
        if weight_col and weight_col in data.columns:
            weights = data[weight_col]
            weight_range = weights.max() - weights.min()
            sizes = (
                40 + (weights - weights.min()) / weight_range * 160
                if weight_range > 0
                else 120
            )
        else:
            sizes = 80
        return cls(
            data=data,
            sizes=sizes,
            study_col=study_col,
            effect_col=effect_col,
            ci_lower_col=ci_lower_col,
            ci_upper_col=ci_upper_col,
            null_value=null_value,
            effect_label=effect_label,
            variant_id=variant_id,
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the confidence intervals, effect markers, null line, and axes."""
        df = self.data
        backend.errorbar_h(
            ax,
            x=df[self.effect_col],
            y=df["y_pos"],
            xerr_lower=df[self.effect_col] - df[self.ci_lower_col],
            xerr_upper=df[self.ci_upper_col] - df[self.effect_col],
            color="black",
            linewidth=1.5,
            capsize=3,
            zorder=2,
        )
        backend.scatter(
            ax,
            df[self.effect_col],
            df["y_pos"],
            colors="#4169E1",
            sizes=self.sizes,
            marker="s",
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )
        backend.axvline(
            ax, x=self.null_value, color="grey", linestyle="--", linewidth=1, alpha=0.7
        )
        backend.set_xlabel(ax, self.effect_label)
        backend.set_ylim(ax, -0.5, len(df) - 0.5)
        x_min = min(df[self.ci_lower_col].min(), self.null_value)
        x_max = max(df[self.ci_upper_col].max(), self.null_value)
        padding = (x_max - x_min) * 0.1
        backend.set_xlim(ax, x_min - padding, x_max + padding)
        backend.set_yticks(
            ax, df["y_pos"].tolist(), df[self.study_col].tolist(), fontsize=10
        )
        backend.set_title(ax, f"Forest Plot: {self.variant_id}")

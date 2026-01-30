"""Statistical visualization plotter for PheWAS and forest plots.

Provides variant-centric visualizations:
- PheWAS plots showing associations across phenotypes
- Forest plots showing effect sizes with confidence intervals
"""

from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd

from ._plotter_utils import transform_pvalues
from .backends import BackendType, get_backend
from .colors import get_phewas_category_palette
from .forest import validate_forest_df
from .phewas import validate_phewas_df

# Default significance threshold
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8


class StatsPlotter:
    """Statistical visualization plotter for PheWAS and forest plots.

    Creates variant-centric visualizations for phenome-wide associations
    and meta-analysis forest plots.

    Args:
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').
        genomewide_threshold: P-value threshold for significance line.

    Example:
        >>> plotter = StatsPlotter()
        >>> fig = plotter.plot_phewas(phewas_df, variant_id="rs12345")
        >>> fig.savefig("phewas.png", dpi=150)
    """

    def __init__(
        self,
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
    ):
        """Initialize the stats plotter."""
        self._backend = get_backend(backend)
        self.genomewide_threshold = genomewide_threshold

    def plot_phewas(
        self,
        phewas_df: pd.DataFrame,
        variant_id: str,
        phenotype_col: str = "phenotype",
        p_col: str = "p_value",
        category_col: str = "category",
        effect_col: Optional[str] = None,
        significance_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
        figsize: Tuple[float, float] = (10, 8),
    ) -> Any:
        """Create a PheWAS (Phenome-Wide Association Study) plot.

        Shows associations of a single variant across multiple phenotypes,
        with phenotypes grouped by category and colored accordingly.

        Args:
            phewas_df: DataFrame with phenotype associations.
            variant_id: Variant identifier (e.g., "rs12345") for plot title.
            phenotype_col: Column name for phenotype names.
            p_col: Column name for p-values.
            category_col: Column name for phenotype categories.
            effect_col: Optional column name for effect direction (beta/OR).
            significance_threshold: P-value threshold for significance line.
            figsize: Figure size as (width, height).

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_phewas(
            ...     phewas_df,
            ...     variant_id="rs12345",
            ...     category_col="category",
            ... )
        """
        validate_phewas_df(phewas_df, phenotype_col, p_col, category_col)

        df = phewas_df.copy()
        df = transform_pvalues(df, p_col)

        # Sort by category then by p-value for consistent ordering
        if category_col in df.columns:
            df = df.sort_values([category_col, p_col])
            categories = df[category_col].unique().tolist()
            palette = get_phewas_category_palette(categories)
        else:
            df = df.sort_values(p_col)
            categories = []
            palette = {}

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Assign y-positions (one per phenotype)
        df["y_pos"] = range(len(df))

        # Plot points by category
        if categories:
            for cat in categories:
                # Handle NaN category: NaN == NaN is False in pandas
                if pd.isna(cat):
                    cat_data = df[df[category_col].isna()]
                else:
                    cat_data = df[df[category_col] == cat]
                # Use upward triangles for positive effects, circles otherwise
                if effect_col and effect_col in cat_data.columns:
                    # Vectorized: split by effect sign, 2 scatter calls per category
                    pos_data = cat_data[cat_data[effect_col] >= 0]
                    neg_data = cat_data[cat_data[effect_col] < 0]

                    if not pos_data.empty:
                        self._backend.scatter(
                            ax,
                            pos_data["neglog10p"],
                            pos_data["y_pos"],
                            colors=palette[cat],
                            sizes=60,
                            marker="^",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
                    if not neg_data.empty:
                        self._backend.scatter(
                            ax,
                            neg_data["neglog10p"],
                            neg_data["y_pos"],
                            colors=palette[cat],
                            sizes=60,
                            marker="v",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
                else:
                    self._backend.scatter(
                        ax,
                        cat_data["neglog10p"],
                        cat_data["y_pos"],
                        colors=palette[cat],
                        sizes=60,
                        marker="o",
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=2,
                    )
        else:
            self._backend.scatter(
                ax,
                df["neglog10p"],
                df["y_pos"],
                colors="#4169E1",
                sizes=60,
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
            )

        # Add significance threshold line
        sig_line = -np.log10(significance_threshold)
        self._backend.axvline(
            ax, x=sig_line, color="red", linestyle="--", linewidth=1, alpha=0.7
        )

        # Set axis labels and limits
        self._backend.set_xlabel(ax, r"$-\log_{10}$ P")
        self._backend.set_ylabel(ax, "Phenotype")
        self._backend.set_ylim(ax, -0.5, len(df) - 0.5)

        # Set y-tick labels to phenotype names
        self._backend.set_yticks(
            ax,
            positions=df["y_pos"].tolist(),
            labels=df[phenotype_col].tolist(),
            fontsize=8,
        )

        self._backend.set_title(ax, f"PheWAS: {variant_id}")
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

    def plot_forest(
        self,
        forest_df: pd.DataFrame,
        variant_id: str,
        study_col: str = "study",
        effect_col: str = "effect",
        ci_lower_col: str = "ci_lower",
        ci_upper_col: str = "ci_upper",
        weight_col: Optional[str] = None,
        null_value: float = 0.0,
        effect_label: str = "Effect Size",
        figsize: Tuple[float, float] = (8, 6),
    ) -> Any:
        """Create a forest plot showing effect sizes with confidence intervals.

        Args:
            forest_df: DataFrame with effect sizes and confidence intervals.
            variant_id: Variant identifier for plot title.
            study_col: Column name for study/phenotype names.
            effect_col: Column name for effect sizes.
            ci_lower_col: Column name for lower confidence interval.
            ci_upper_col: Column name for upper confidence interval.
            weight_col: Optional column for study weights (affects marker size).
            null_value: Reference value for null effect (0 for beta, 1 for OR).
            effect_label: X-axis label.
            figsize: Figure size as (width, height).

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_forest(
            ...     forest_df,
            ...     variant_id="rs12345",
            ...     effect_label="Odds Ratio",
            ...     null_value=1.0,
            ... )
        """
        validate_forest_df(forest_df, study_col, effect_col, ci_lower_col, ci_upper_col)

        df = forest_df.copy()

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Assign y-positions (reverse so first study is at top)
        df["y_pos"] = range(len(df) - 1, -1, -1)

        # Calculate marker sizes from weights
        if weight_col and weight_col in df.columns:
            # Scale weights to marker sizes (min 40, max 200)
            weights = df[weight_col]
            min_size, max_size = 40, 200
            weight_range = weights.max() - weights.min()
            if weight_range > 0:
                sizes = min_size + (weights - weights.min()) / weight_range * (
                    max_size - min_size
                )
            else:
                sizes = (min_size + max_size) / 2
        else:
            sizes = 80

        # Calculate error bar extents
        xerr_lower = df[effect_col] - df[ci_lower_col]
        xerr_upper = df[ci_upper_col] - df[effect_col]

        # Plot error bars (confidence intervals)
        self._backend.errorbar_h(
            ax,
            x=df[effect_col],
            y=df["y_pos"],
            xerr_lower=xerr_lower,
            xerr_upper=xerr_upper,
            color="black",
            linewidth=1.5,
            capsize=3,
            zorder=2,
        )

        # Plot effect size markers
        self._backend.scatter(
            ax,
            df[effect_col],
            df["y_pos"],
            colors="#4169E1",
            sizes=sizes,
            marker="s",  # square markers typical for forest plots
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )

        # Add null effect line
        self._backend.axvline(
            ax, x=null_value, color="grey", linestyle="--", linewidth=1, alpha=0.7
        )

        # Set axis labels and limits
        self._backend.set_xlabel(ax, effect_label)
        self._backend.set_ylim(ax, -0.5, len(df) - 0.5)

        # Ensure x-axis includes the null value with some padding
        x_min = min(df[ci_lower_col].min(), null_value)
        x_max = max(df[ci_upper_col].max(), null_value)
        x_padding = (x_max - x_min) * 0.1
        self._backend.set_xlim(ax, x_min - x_padding, x_max + x_padding)

        # Set y-tick labels to study names
        self._backend.set_yticks(
            ax,
            positions=df["y_pos"].tolist(),
            labels=df[study_col].tolist(),
            fontsize=10,
        )

        self._backend.set_title(ax, f"Forest Plot: {variant_id}")
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

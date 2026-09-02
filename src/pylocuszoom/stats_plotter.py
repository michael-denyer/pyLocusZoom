"""Statistical visualization plotter for PheWAS and forest plots.

Provides variant-centric visualizations:
- PheWAS plots showing associations across phenotypes
- Forest plots showing effect sizes with confidence intervals
"""

from typing import Any, Optional, Tuple

import pandas as pd

from ._data import prepare_pvalue_data
from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from ._stats_renderer import StatsRenderer
from .backends import BackendType, get_backend
from .schemas import validate_forest_df, validate_phewas_df


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
        self._renderer = StatsRenderer(self._backend)
        self.genomewide_threshold = genomewide_threshold

    def plot_phewas(
        self,
        phewas_df: pd.DataFrame,
        variant_id: str,
        phenotype_col: str = "phenotype",
        p_col: str = "p_value",
        category_col: str = "category",
        effect_col: Optional[str] = None,
        significance_threshold: ThresholdArg = UNSET,
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
            significance_threshold: P-value threshold for the significance line.
                Defaults to the plotter's ``genomewide_threshold``; pass None to
                draw no line.
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
        significance_threshold = resolve_threshold(
            significance_threshold, self.genomewide_threshold
        )
        validate_phewas_df(phewas_df, phenotype_col, p_col, category_col)

        df = phewas_df.copy()
        df = prepare_pvalue_data(df, p_col)

        return self._renderer.render_phewas(
            df,
            phenotype_col=phenotype_col,
            p_col=p_col,
            category_col=category_col,
            variant_id=variant_id,
            effect_col=effect_col,
            significance_threshold=significance_threshold,
            figsize=figsize,
        )

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

        return self._renderer.render_forest(
            df,
            study_col=study_col,
            effect_col=effect_col,
            ci_lower_col=ci_lower_col,
            ci_upper_col=ci_upper_col,
            weight_col=weight_col,
            null_value=null_value,
            effect_label=effect_label,
            variant_id=variant_id,
            figsize=figsize,
        )

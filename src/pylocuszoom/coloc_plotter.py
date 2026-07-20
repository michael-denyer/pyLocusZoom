"""Colocalization scatter plot for GWAS-eQTL visualization.

Creates scatter plots comparing GWAS -log10(p) vs eQTL -log10(p)
with points colored by LD to the lead SNP.
"""

from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd

from ._coloc_renderer import ColocRenderer
from ._data import P_VALUE_FLOOR
from .backends import BackendType, get_backend
from .coloc import validate_coloc_eqtl_df, validate_coloc_gwas_df
from .colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    LD_NA_COLOR,
    get_ld_color,
)
from .config import ColocConfig


def _resolve_merged_column(
    merged: pd.DataFrame,
    col: Optional[str],
    suffix: str,
) -> Optional[str]:
    """Resolve column name after DataFrame merge.

    When merging DataFrames, pandas adds suffixes to duplicate columns.
    This helper finds the actual column name in the merged DataFrame.

    Args:
        merged: Merged DataFrame to search.
        col: Original column name (or None).
        suffix: Suffix added by merge (e.g., "_gwas" or "_eqtl").

    Returns:
        Resolved column name, or None if col was None or not found.
    """
    if col is None:
        return None
    suffixed = f"{col}{suffix}"
    if suffixed in merged.columns:
        return suffixed
    if col in merged.columns:
        return col
    return None


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


class ColocPlotter:
    """Colocalization scatter plot generator.

    Creates scatter plots comparing GWAS -log10(p) vs eQTL -log10(p)
    with points colored by LD to the lead SNP.

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').

    Example:
        >>> plotter = ColocPlotter()
        >>> fig = plotter.plot_coloc(gwas_df, eqtl_df, lead_snp="rs12345")
        >>> fig.savefig("coloc.png", dpi=150)
    """

    def __init__(
        self,
        backend: BackendType = "matplotlib",
    ):
        """Initialize the colocalization plotter."""
        self._backend = get_backend(backend)
        self._renderer = ColocRenderer(self._backend)
        self.backend_name = backend

    def plot_coloc(
        self,
        gwas_df: pd.DataFrame,
        eqtl_df: pd.DataFrame,
        pos_col: str = "pos",
        gwas_p_col: str = "p_gwas",
        eqtl_p_col: str = "p_eqtl",
        rs_col: Optional[str] = "rs",
        ld_col: Optional[str] = None,
        lead_snp: Optional[str] = None,
        gwas_threshold: float = 5e-8,
        eqtl_threshold: float = 1e-5,
        show_correlation: bool = True,
        color_by_effect: bool = False,
        gwas_effect_col: Optional[str] = None,
        eqtl_effect_col: Optional[str] = None,
        h4_posterior: Optional[float] = None,
        figsize: Tuple[float, float] = (8.0, 8.0),
        title: Optional[str] = None,
    ) -> Any:
        """Create GWAS-eQTL colocalization scatter plot.

        Args:
            gwas_df: GWAS results DataFrame with positions and p-values.
            eqtl_df: eQTL results DataFrame with positions and p-values.
            pos_col: Column name for genomic positions (must exist in both).
            gwas_p_col: Column name for GWAS p-values.
            eqtl_p_col: Column name for eQTL p-values.
            rs_col: Column name for SNP IDs (optional, for labeling lead SNP).
            ld_col: Column name for LD R² values in GWAS df (optional).
            lead_snp: SNP ID to highlight as lead variant. If None and ld_col
                is provided, auto-selects SNP with highest combined -log10(p).
            gwas_threshold: P-value threshold for GWAS significance line.
            eqtl_threshold: P-value threshold for eQTL significance line.
            show_correlation: Whether to display Pearson correlation.
            color_by_effect: Whether to color points by effect direction agreement.
            gwas_effect_col: Column name for GWAS effect sizes (required if
                color_by_effect=True).
            eqtl_effect_col: Column name for eQTL effect sizes (required if
                color_by_effect=True).
            h4_posterior: Optional COLOC H4 posterior probability to display.
            figsize: Figure size as (width, height).
            title: Plot title.

        Returns:
            Figure object (type depends on backend).

        Raises:
            ValidationError: If required columns are missing or invalid.
            ValueError: If no overlapping positions between GWAS and eQTL.
            ValueError: If lead_snp specified but not found in merged data.
            ValueError: If color_by_effect=True but effect columns not provided.
            ValueError: If h4_posterior is not in [0, 1] range.

        Example:
            >>> fig = plotter.plot_coloc(
            ...     gwas_df, eqtl_df,
            ...     ld_col="ld", lead_snp="rs12345",
            ... )
            >>> # With effect coloring
            >>> fig = plotter.plot_coloc(
            ...     gwas_df, eqtl_df,
            ...     color_by_effect=True,
            ...     gwas_effect_col="beta_gwas",
            ...     eqtl_effect_col="beta_eqtl",
            ... )
        """
        config = ColocConfig(
            pos_col=pos_col,
            gwas_p_col=gwas_p_col,
            eqtl_p_col=eqtl_p_col,
            rs_col=rs_col,
            ld_col=ld_col,
            lead_snp=lead_snp,
            gwas_threshold=gwas_threshold,
            eqtl_threshold=eqtl_threshold,
            show_correlation=show_correlation,
            color_by_effect=color_by_effect,
            gwas_effect_col=gwas_effect_col,
            eqtl_effect_col=eqtl_effect_col,
            h4_posterior=h4_posterior,
            figsize=figsize,
        )
        # Validate inputs
        validate_coloc_gwas_df(
            gwas_df,
            config.pos_col,
            config.gwas_p_col,
            config.rs_col,
        )
        validate_coloc_eqtl_df(
            eqtl_df,
            config.pos_col,
            config.eqtl_p_col,
            config.rs_col,
        )

        # Merge DataFrames on position
        merged = pd.merge(
            gwas_df,
            eqtl_df,
            on=config.pos_col,
            how="inner",
            suffixes=("_gwas", "_eqtl"),
        )

        if len(merged) == 0:
            raise ValueError(
                "No overlapping positions between GWAS and eQTL DataFrames"
            )

        # Resolve column names after merge (pandas adds suffixes to duplicates)
        merged_rs_col = _resolve_merged_column(merged, config.rs_col, "_gwas")
        ld_col_merged = _resolve_merged_column(merged, config.ld_col, "_gwas")
        gwas_p_merged = _resolve_merged_column(merged, config.gwas_p_col, "_gwas")
        eqtl_p_merged = _resolve_merged_column(merged, config.eqtl_p_col, "_eqtl")

        # Coloc transforms two merged p-value columns at once, so it does its
        # own -log10 here rather than the single-column prepare_pvalue_data
        # intake; the p-value floor stays shared via P_VALUE_FLOOR.
        merged["neglog10_gwas"] = -np.log10(
            merged[gwas_p_merged].clip(lower=P_VALUE_FLOOR)
        )
        merged["neglog10_eqtl"] = -np.log10(
            merged[eqtl_p_merged].clip(lower=P_VALUE_FLOOR)
        )

        # Drop rows with NaN in either transformed p-value
        merged = merged.dropna(subset=["neglog10_gwas", "neglog10_eqtl"])

        if len(merged) == 0:
            raise ValueError("No valid data points after removing NaN p-values")

        # Resolve effect columns if coloring by effect direction
        gwas_effect_merged = None
        eqtl_effect_merged = None
        if config.color_by_effect:
            gwas_effect_merged = _resolve_merged_column(
                merged, config.gwas_effect_col, "_gwas"
            )
            if gwas_effect_merged is None:
                raise ValueError(
                    "gwas_effect_col "
                    f"'{config.gwas_effect_col}' not found in merged data"
                )
            eqtl_effect_merged = _resolve_merged_column(
                merged, config.eqtl_effect_col, "_eqtl"
            )
            if eqtl_effect_merged is None:
                raise ValueError(
                    "eqtl_effect_col "
                    f"'{config.eqtl_effect_col}' not found in merged data"
                )

        # Apply coloring based on mode
        if config.color_by_effect:
            merged["color"] = merged.apply(
                lambda row: _get_effect_agreement_color(
                    row[gwas_effect_merged], row[eqtl_effect_merged]
                ),
                axis=1,
            )
        elif ld_col_merged is not None:
            merged["color"] = merged[ld_col_merged].apply(get_ld_color)
        else:
            merged["color"] = LD_NA_COLOR

        # Determine lead SNP index
        lead_idx = None
        if config.lead_snp is not None:
            if merged_rs_col is None:
                raise ValueError(
                    f"lead_snp '{config.lead_snp}' specified but rs_col not found"
                )
            matches = merged[merged[merged_rs_col] == config.lead_snp]
            if len(matches) == 0:
                raise ValueError(
                    f"lead_snp '{config.lead_snp}' not found in merged data"
                )
            lead_idx = matches.index[0]
        elif ld_col_merged is not None:
            # Auto-select: highest combined -log10(p)
            merged["combined_score"] = merged["neglog10_gwas"] + merged["neglog10_eqtl"]
            lead_idx = merged["combined_score"].idxmax()

        return self._renderer.render(
            merged,
            lead_idx=lead_idx,
            merged_rs_col=merged_rs_col,
            ld_col_merged=ld_col_merged,
            gwas_threshold=config.gwas_threshold,
            eqtl_threshold=config.eqtl_threshold,
            show_correlation=config.show_correlation,
            color_by_effect=config.color_by_effect,
            h4_posterior=config.h4_posterior,
            title=title,
            figsize=config.figsize,
        )

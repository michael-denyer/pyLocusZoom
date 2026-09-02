"""Colocalization scatter plot for GWAS-eQTL visualization.

Creates scatter plots comparing GWAS -log10(p) vs eQTL -log10(p)
with points colored by LD to the lead SNP.
"""

from dataclasses import dataclass
from typing import Any, Optional, Tuple

import pandas as pd

from ._coloc_renderer import ColocRequest, render_coloc
from ._data import prepare_pvalue_data
from ._plotter_utils import (
    DEFAULT_EQTL_THRESHOLD,
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from .backends import BackendType, get_backend
from .colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    LD_NA_COLOR,
    get_ld_color,
)
from .config import ColocConfig
from .schemas import validate_coloc_df


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


@dataclass(frozen=True)
class _MergedColoc:
    """A GWAS/eQTL frame merged on position, with column names resolved.

    The merge suffixes any column name the two frames share, so every later
    step needs the post-merge name rather than the caller's. Carrying them
    beside the frame keeps that resolution in one place.

    Attributes:
        data: Merged rows carrying ``neglog10_gwas`` and ``neglog10_eqtl``.
        rs_col: Post-merge SNP ID column, or None if the frames had none.
        ld_col: Post-merge LD column, or None if the caller supplied none.
    """

    data: pd.DataFrame
    rs_col: Optional[str]
    ld_col: Optional[str]


def _merge_and_transform(
    gwas_df: pd.DataFrame,
    eqtl_df: pd.DataFrame,
    config: ColocConfig,
) -> _MergedColoc:
    """Transform both p-value columns through the shared intake, then merge.

    Args:
        gwas_df: Validated GWAS results.
        eqtl_df: Validated eQTL results.
        config: Validated plot configuration.

    Returns:
        The merged frame and its resolved column names.

    Raises:
        ValueError: If the frames share no positions.
    """
    merged = pd.merge(
        prepare_pvalue_data(gwas_df, config.gwas_p_col, out_col="neglog10_gwas"),
        prepare_pvalue_data(eqtl_df, config.eqtl_p_col, out_col="neglog10_eqtl"),
        on=config.pos_col,
        how="inner",
        suffixes=("_gwas", "_eqtl"),
    )
    if len(merged) == 0:
        raise ValueError("No overlapping positions between GWAS and eQTL DataFrames")

    return _MergedColoc(
        data=merged,
        rs_col=_resolve_merged_column(merged, config.rs_col, "_gwas"),
        ld_col=_resolve_merged_column(merged, config.ld_col, "_gwas"),
    )


def _assign_colors(merged: _MergedColoc, config: ColocConfig) -> pd.Series:
    """Colour every merged point by effect agreement, by LD, or not at all.

    Args:
        merged: Output of ``_merge_and_transform``.
        config: Validated plot configuration.

    Returns:
        One hex colour per row, aligned to ``merged.data``.

    Raises:
        ValueError: If ``color_by_effect`` is set but either effect column is
            absent from the merged frame.
    """
    if config.color_by_effect:
        gwas_effect = _resolve_merged_column(
            merged.data, config.gwas_effect_col, "_gwas"
        )
        if gwas_effect is None:
            raise ValueError(
                f"gwas_effect_col '{config.gwas_effect_col}' not found in merged data"
            )
        eqtl_effect = _resolve_merged_column(
            merged.data, config.eqtl_effect_col, "_eqtl"
        )
        if eqtl_effect is None:
            raise ValueError(
                f"eqtl_effect_col '{config.eqtl_effect_col}' not found in merged data"
            )
        return merged.data.apply(
            lambda row: _get_effect_agreement_color(row[gwas_effect], row[eqtl_effect]),
            axis=1,
        )
    if merged.ld_col is not None:
        return merged.data[merged.ld_col].apply(get_ld_color)
    return pd.Series(LD_NA_COLOR, index=merged.data.index)


def _resolve_lead_idx(merged: _MergedColoc, config: ColocConfig) -> Optional[Any]:
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
        ValueError: If ``lead_snp`` is named but the merged frame has no SNP
            ID column or no row matching it.
    """
    if config.lead_snp is not None:
        if merged.rs_col is None:
            raise ValueError(
                f"lead_snp '{config.lead_snp}' specified but rs_col not found"
            )
        matches = merged.data[merged.data[merged.rs_col] == config.lead_snp]
        if len(matches) == 0:
            raise ValueError(f"lead_snp '{config.lead_snp}' not found in merged data")
        return matches.index[0]
    if merged.ld_col is not None:
        combined = merged.data["neglog10_gwas"] + merged.data["neglog10_eqtl"]
        return combined.idxmax()
    return None


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
        genomewide_threshold: P-value threshold for the GWAS significance line.
        eqtl_threshold: P-value threshold for the eQTL significance line.

    Example:
        >>> plotter = ColocPlotter()
        >>> fig = plotter.plot_coloc(gwas_df, eqtl_df, lead_snp="rs12345")
        >>> fig.savefig("coloc.png", dpi=150)
    """

    def __init__(
        self,
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
        eqtl_threshold: float = DEFAULT_EQTL_THRESHOLD,
    ):
        """Initialize the colocalization plotter."""
        self._backend = get_backend(backend)
        self.genomewide_threshold = genomewide_threshold
        self.eqtl_threshold = eqtl_threshold

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
        gwas_threshold: ThresholdArg = UNSET,
        eqtl_threshold: ThresholdArg = UNSET,
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
            gwas_threshold: Significance threshold for the GWAS line. Defaults
                to the plotter's ``genomewide_threshold``; pass None to draw no
                line.
            eqtl_threshold: Significance threshold for the eQTL line. Defaults
                to the plotter's ``eqtl_threshold``; pass None to draw no line.
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
            gwas_threshold=resolve_threshold(gwas_threshold, self.genomewide_threshold),
            eqtl_threshold=resolve_threshold(eqtl_threshold, self.eqtl_threshold),
            show_correlation=show_correlation,
            color_by_effect=color_by_effect,
            gwas_effect_col=gwas_effect_col,
            eqtl_effect_col=eqtl_effect_col,
            h4_posterior=h4_posterior,
            figsize=figsize,
        )
        validate_coloc_df(
            gwas_df,
            "GWAS DataFrame",
            config.pos_col,
            config.gwas_p_col,
            config.rs_col,
        )
        validate_coloc_df(
            eqtl_df,
            "eQTL DataFrame",
            config.pos_col,
            config.eqtl_p_col,
            config.rs_col,
        )

        merged = _merge_and_transform(gwas_df, eqtl_df, config)
        merged.data["color"] = _assign_colors(merged, config)
        lead_idx = _resolve_lead_idx(merged, config)

        return render_coloc(
            self._backend,
            ColocRequest(
                merged=merged.data,
                config=config,
                rs_col=merged.rs_col,
                ld_col=merged.ld_col,
                lead_idx=lead_idx,
                title=title,
            ),
        )

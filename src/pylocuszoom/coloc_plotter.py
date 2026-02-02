"""Colocalization scatter plot for GWAS-eQTL visualization.

Creates scatter plots comparing GWAS -log10(p) vs eQTL -log10(p)
with points colored by LD to the lead SNP.
"""

from typing import Any, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from .backends import BackendType, get_backend
from .coloc import validate_coloc_eqtl_df, validate_coloc_gwas_df
from .colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    LD_BINS,
    LD_NA_COLOR,
    LEAD_SNP_COLOR,
    get_ld_color,
)


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
        # Validate inputs
        validate_coloc_gwas_df(gwas_df, pos_col, gwas_p_col, rs_col)
        validate_coloc_eqtl_df(eqtl_df, pos_col, eqtl_p_col, rs_col)

        # Validate effect coloring parameters
        if color_by_effect:
            if gwas_effect_col is None or eqtl_effect_col is None:
                raise ValueError(
                    "color_by_effect=True requires gwas_effect_col and eqtl_effect_col"
                )

        # Validate h4_posterior range
        if h4_posterior is not None:
            if not (0 <= h4_posterior <= 1):
                raise ValueError(f"h4_posterior must be in [0, 1], got {h4_posterior}")

        # Merge DataFrames on position
        merged = pd.merge(
            gwas_df,
            eqtl_df,
            on=pos_col,
            how="inner",
            suffixes=("_gwas", "_eqtl"),
        )

        if len(merged) == 0:
            raise ValueError(
                "No overlapping positions between GWAS and eQTL DataFrames"
            )

        # Resolve column names after merge (pandas adds suffixes to duplicates)
        merged_rs_col = _resolve_merged_column(merged, rs_col, "_gwas")
        ld_col_merged = _resolve_merged_column(merged, ld_col, "_gwas")
        gwas_p_merged = _resolve_merged_column(merged, gwas_p_col, "_gwas")
        eqtl_p_merged = _resolve_merged_column(merged, eqtl_p_col, "_eqtl")

        # Transform p-values to -log10(p)
        merged["neglog10_gwas"] = -np.log10(merged[gwas_p_merged].clip(lower=1e-300))
        merged["neglog10_eqtl"] = -np.log10(merged[eqtl_p_merged].clip(lower=1e-300))

        # Drop rows with NaN in either transformed p-value
        merged = merged.dropna(subset=["neglog10_gwas", "neglog10_eqtl"])

        if len(merged) == 0:
            raise ValueError("No valid data points after removing NaN p-values")

        # Resolve effect columns if coloring by effect direction
        gwas_effect_merged = None
        eqtl_effect_merged = None
        if color_by_effect:
            gwas_effect_merged = _resolve_merged_column(
                merged, gwas_effect_col, "_gwas"
            )
            if gwas_effect_merged is None:
                raise ValueError(
                    f"gwas_effect_col '{gwas_effect_col}' not found in merged data"
                )
            eqtl_effect_merged = _resolve_merged_column(
                merged, eqtl_effect_col, "_eqtl"
            )
            if eqtl_effect_merged is None:
                raise ValueError(
                    f"eqtl_effect_col '{eqtl_effect_col}' not found in merged data"
                )

        # Apply coloring based on mode
        if color_by_effect:
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
        if lead_snp is not None:
            if merged_rs_col is None:
                raise ValueError(
                    f"lead_snp '{lead_snp}' specified but rs_col not found"
                )
            matches = merged[merged[merged_rs_col] == lead_snp]
            if len(matches) == 0:
                raise ValueError(f"lead_snp '{lead_snp}' not found in merged data")
            lead_idx = matches.index[0]
        elif ld_col_merged is not None:
            # Auto-select: highest combined -log10(p)
            merged["combined_score"] = merged["neglog10_gwas"] + merged["neglog10_eqtl"]
            lead_idx = merged["combined_score"].idxmax()

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Separate lead SNP from other points
        if lead_idx is not None:
            lead_row = merged.loc[[lead_idx]]
            other_rows = merged.drop(lead_idx)
        else:
            lead_row = pd.DataFrame()
            other_rows = merged

        # Plot non-lead points
        if len(other_rows) > 0:
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

        # Plot lead SNP as diamond
        if len(lead_row) > 0:
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

            # Add lead SNP label
            if merged_rs_col is not None:
                label = lead_row[merged_rs_col].values[0]
                x_pos = lead_row["neglog10_gwas"].values[0]
                y_pos = lead_row["neglog10_eqtl"].values[0]
                self._backend.add_text(
                    ax,
                    x_pos,
                    y_pos + 0.5,  # Offset above the point
                    label,
                    fontsize=9,
                    ha="center",
                    va="bottom",
                )

        # Add significance threshold lines
        gwas_sig_line = -np.log10(gwas_threshold)
        eqtl_sig_line = -np.log10(eqtl_threshold)

        self._backend.axvline(
            ax, x=gwas_sig_line, color="grey", linestyle="--", linewidth=1, alpha=0.7
        )
        self._backend.axhline(
            ax, y=eqtl_sig_line, color="grey", linestyle="--", linewidth=1, alpha=0.7
        )

        # Calculate data bounds once for text positioning
        x_min, x_max = merged["neglog10_gwas"].min(), merged["neglog10_gwas"].max()
        y_min, y_max = merged["neglog10_eqtl"].min(), merged["neglog10_eqtl"].max()
        x_range = x_max - x_min
        y_range = y_max - y_min

        # Display correlation in top-left corner
        if show_correlation and len(merged) >= 3:
            r, p = stats.pearsonr(merged["neglog10_gwas"], merged["neglog10_eqtl"])
            p_str = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
            corr_text = f"r = {r:.3f}\n{p_str}"
            self._backend.add_text(
                ax,
                x_min + 0.05 * x_range,
                y_max - 0.05 * y_range,
                corr_text,
                fontsize=10,
                ha="left",
                va="top",
            )

        # Display H4 posterior probability in bottom-right corner
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

        # Set axis labels
        self._backend.set_xlabel(ax, r"GWAS $-\log_{10}$ P")
        self._backend.set_ylabel(ax, r"eQTL $-\log_{10}$ P")

        # Set title
        if title:
            self._backend.set_title(ax, title)

        # Hide top and right spines
        self._backend.hide_spines(ax, ["top", "right"])

        # Add legend
        if color_by_effect:
            self._add_effect_legend(ax)
        elif ld_col_merged is not None:
            self._backend.add_ld_legend(ax, LD_BINS, LEAD_SNP_COLOR)

        # Finalize layout
        self._backend.finalize_layout(fig)

        return fig

    def _add_effect_legend(self, ax: Any) -> None:
        """Add effect direction legend to plot (all backends)."""
        effect_bins = [
            (0.0, "Same direction", EFFECT_CONGRUENT_COLOR),
            (0.0, "Opposite direction", EFFECT_INCONGRUENT_COLOR),
            (0.0, "Missing effect", LD_NA_COLOR),
        ]
        self._backend.add_effect_legend(ax, effect_bins)

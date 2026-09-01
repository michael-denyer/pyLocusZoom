"""SNP label placement for regional association plots.

Provides automatic labeling of top significant SNPs with:
- SNP ID (rs number)
- Automatic overlap avoidance (if adjustText installed)
"""

from typing import List, Optional

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.text import Annotation

from pylocuszoom.logging import logger


def add_snp_labels(
    ax: Axes,
    df: pd.DataFrame,
    pos_col: str = "ps",
    neglog10p_col: str = "neglog10p",
    rs_col: str = "rs",
    label_top_n: int = 5,
    max_label_length: int = 15,
    adjust: bool = True,
    lead_pos: Optional[int] = None,
    region_span: Optional[int] = None,
    min_label_distance: float = 0.05,
) -> List[Annotation]:
    """Add text labels to top SNPs in the regional plot.

    Labels the most significant SNPs with their SNP ID (rs number).
    When a lead SNP position is provided, nearby non-lead SNPs are
    excluded to avoid overlapping connector lines.

    Args:
        ax: Matplotlib axes object.
        df: DataFrame with SNP data. Must have the specified position,
            neglog10p, and rs columns.
        pos_col: Column name for position.
        neglog10p_col: Column name for -log10(p-value).
        rs_col: Column name for SNP ID.
        label_top_n: Number of top SNPs to label.
        max_label_length: Maximum label length before truncation.
        adjust: If True, run adjustText immediately. If False, caller must
            call adjust_snp_labels() after setting axis limits.
        lead_pos: Position of the lead/index SNP. Non-lead SNPs closer
            than ``min_label_distance * region_span`` base pairs are skipped.
            The lead SNP is identified by exact position match.
        region_span: Width of the visible region in base pairs. Must be
            positive for filtering to take effect. If not provided when
            ``lead_pos`` is set, a warning is logged and filtering is
            skipped.
        min_label_distance: Minimum distance from lead SNP as a fraction
            of ``region_span``. Non-lead SNPs closer than this are not
            labeled. Defaults to 0.05 (5%). Must be between 0 and 1.

    Returns:
        List of matplotlib text annotation objects.

    Example:
        >>> fig, ax = plt.subplots()
        >>> # ... plot your data ...
        >>> texts = add_snp_labels(ax, df, label_top_n=5)
    """
    if neglog10p_col not in df.columns:
        raise ValueError(
            f"Column '{neglog10p_col}' not found in DataFrame. "
            "Ensure -log10(p) values are calculated before calling add_snp_labels."
        )

    # Filter eligible SNPs FIRST (drop near-lead non-lead rows), then take top N.
    # Otherwise nlargest can be saturated by the lead's neighbors and the mask
    # leaves only the lead itself, producing far fewer labels than requested.
    eligible = df
    if lead_pos is not None and region_span is not None and region_span > 0:
        if not 0 <= min_label_distance <= 1:
            raise ValueError(
                f"min_label_distance must be between 0 and 1, got {min_label_distance}"
            )
        min_dist_bp = min_label_distance * region_span
        mask = (df[pos_col] == lead_pos) | (
            (df[pos_col] - lead_pos).abs() >= min_dist_bp
        )
        eligible = df[mask]
    elif lead_pos is not None and (region_span is None or region_span <= 0):
        logger.warning(
            "lead_pos provided without valid region_span — "
            "proximity filtering disabled. Pass region_span=(end - start) to enable."
        )

    top_snps = eligible.nlargest(label_top_n, neglog10p_col)

    texts = []
    used_labels = set()  # Track used labels to avoid duplicates

    for _, snp in top_snps.iterrows():
        x = snp[pos_col]
        y = snp[neglog10p_col]

        # Use SNP ID as label
        label = str(snp[rs_col])

        # Skip duplicate labels
        if label in used_labels:
            continue
        used_labels.add(label)

        # Truncate long labels
        if len(label) > max_label_length:
            label = label[: max_label_length - 3] + "..."

        # Add text annotation centered above marker
        text = ax.annotate(
            label,
            xy=(x, y),
            xytext=(0, 7),
            textcoords="offset points",
            fontsize=6,
            fontweight="bold",
            color="#333333",
            ha="center",
            va="bottom",
            zorder=15,
            bbox=dict(
                boxstyle="round,pad=0.2",
                facecolor="white",
                edgecolor="none",
                alpha=0.8,
            ),
        )
        texts.append(text)

    if adjust:
        adjust_snp_labels(ax, texts)

    return texts


def adjust_snp_labels(ax: Axes, texts: List[Annotation]) -> None:
    """Adjust SNP label positions to avoid overlaps.

    This function should be called AFTER all axis limits have been set,
    as adjustText needs to know the final plot bounds to position labels
    correctly within the visible area.

    Args:
        ax: Matplotlib axes object.
        texts: List of text annotation objects from add_snp_labels().

    Example:
        >>> texts = add_snp_labels(ax, df, adjust=False)
        >>> ax.set_xlim(start, end)
        >>> ax.set_ylim(0, max_y)
        >>> adjust_snp_labels(ax, texts)
    """
    if len(texts) <= 1:
        return

    try:
        from adjustText import adjust_text

        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="none", lw=0),
            expand_points=(1.5, 1.5),
        )
    except ImportError:
        logger.warning(
            "adjustText not installed - SNP labels may overlap. "
            "Install with: pip install adjustText"
        )

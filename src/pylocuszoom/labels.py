"""SNP label placement for regional association plots.

Provides automatic labeling of top significant SNPs with:
- SNP ID (rs number)
- Automatic overlap avoidance (if adjustText installed)
"""

from typing import Any, List, Optional, Union

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.text import Annotation


def add_snp_labels(
    ax: Axes,
    df: pd.DataFrame,
    pos_col: str = "ps",
    neglog10p_col: str = "neglog10p",
    rs_col: str = "rs",
    label_top_n: int = 5,
    genes_df: Optional[pd.DataFrame] = None,
    chrom: Optional[Union[int, str]] = None,
    max_label_length: int = 15,
    **kwargs: Any,
) -> List[Annotation]:
    """Add text labels to top SNPs in the regional plot.

    Labels the most significant SNPs with their SNP ID (rs number).

    Args:
        ax: Matplotlib axes object.
        df: DataFrame with SNP data. Must have the specified position,
            neglog10p, and rs columns.
        pos_col: Column name for position.
        neglog10p_col: Column name for -log10(p-value).
        rs_col: Column name for SNP ID.
        label_top_n: Number of top SNPs to label.
        genes_df: Unused, kept for backward compatibility.
        chrom: Unused, kept for backward compatibility.
        max_label_length: Maximum label length before truncation.

    Returns:
        List of matplotlib text annotation objects.

    Example:
        >>> fig, ax = plt.subplots()
        >>> # ... plot your data ...
        >>> texts = add_snp_labels(ax, df, label_top_n=5)
    """
    # genes_df and chrom are unused but kept for backward compatibility
    del genes_df, chrom, kwargs
    if neglog10p_col not in df.columns:
        raise ValueError(
            f"Column '{neglog10p_col}' not found in DataFrame. "
            "Ensure -log10(p) values are calculated before calling add_snp_labels."
        )

    # Get top N SNPs by -log10(p)
    top_snps = df.nlargest(label_top_n, neglog10p_col)

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

    # Only use adjustText when there are multiple labels to avoid overlap
    if len(texts) > 1:
        try:
            from adjustText import adjust_text

            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                expand_points=(1.5, 1.5),
            )
        except ImportError:
            # adjustText not installed, labels may overlap
            pass

    return texts

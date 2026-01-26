"""SNP label placement for regional association plots.

Provides automatic labeling of top significant SNPs with:
- SNP ID (rs number)
- Nearest gene name (if gene annotations provided)
- Automatic overlap avoidance (if adjustText installed)
"""

from typing import List, Optional, Union

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.text import Annotation

from .gene_track import get_nearest_gene


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
) -> List[Annotation]:
    """Add text labels to top SNPs in the regional plot.

    Labels the most significant SNPs with either their SNP ID
    or the nearest gene name (if genes_df provided).

    Args:
        ax: Matplotlib axes object.
        df: DataFrame with SNP data. Must have the specified position,
            neglog10p, and rs columns.
        pos_col: Column name for position.
        neglog10p_col: Column name for -log10(p-value).
        rs_col: Column name for SNP ID.
        label_top_n: Number of top SNPs to label.
        genes_df: Optional gene annotations for gene-based labels.
            If provided with chrom, labels will show nearest gene name
            instead of SNP ID.
        chrom: Chromosome number. Required if genes_df is provided.
        max_label_length: Maximum label length before truncation.

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

    # Get top N SNPs by -log10(p)
    top_snps = df.nlargest(label_top_n, neglog10p_col)

    texts = []
    for _, snp in top_snps.iterrows():
        x = snp[pos_col]
        y = snp[neglog10p_col]

        # Determine label text
        label = str(snp[rs_col])

        # Try to get gene name if genes_df provided
        if genes_df is not None and chrom is not None:
            nearest_gene = get_nearest_gene(genes_df, chrom, int(x))
            if nearest_gene:
                label = nearest_gene

        # Truncate long labels
        if len(label) > max_label_length:
            label = label[: max_label_length - 3] + "..."

        # Add text annotation with offset
        text = ax.annotate(
            label,
            xy=(x, y),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            color="#333333",
            ha="left",
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

    # Try to adjust text positions to avoid overlap
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

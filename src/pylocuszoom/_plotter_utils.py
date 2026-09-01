"""Shared utilities for plotter classes.

Internal module - not part of public API.
"""

from typing import Any, Optional, Union

import numpy as np
import pandas as pd

# Significance thresholds
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8


class _Unset:
    """Type of the ``UNSET`` sentinel."""

    def __repr__(self) -> str:
        return "UNSET"


# Marks a threshold argument the caller did not supply, so it falls back to the
# plotter's genomewide_threshold. None cannot carry that meaning: it already
# means "draw no significance line", and a caller must keep being able to say so.
UNSET = _Unset()

# A significance threshold as a caller may express it: a p-value, None for no
# line, or UNSET to inherit the plotter's own threshold.
ThresholdArg = Union[float, None, _Unset]


def resolve_threshold(
    supplied: ThresholdArg, plotter_default: Optional[float]
) -> Optional[float]:
    """Pick the threshold to draw at, honouring an explicit None.

    Args:
        supplied: The caller's threshold argument.
        plotter_default: The plotter's own ``genomewide_threshold``.

    Returns:
        The p-value to draw the line at, or None to draw no line.
    """
    return plotter_default if isinstance(supplied, _Unset) else supplied


# Manhattan/QQ plot styling constants
MANHATTAN_POINT_SIZE = 10
MANHATTAN_CATEGORICAL_POINT_SIZE = 30
QQ_POINT_SIZE = 10
POINT_EDGE_COLOR = "black"
MANHATTAN_EDGE_WIDTH = 0.1
QQ_EDGE_WIDTH = 0.02
QQ_POINT_COLOR = "#1f77b4"
QQ_CI_COLOR = "#CCCCCC"
QQ_CI_ALPHA = 0.5
SIGNIFICANCE_LINE_COLOR = "red"


def add_significance_line(
    backend: Any,
    ax: Any,
    threshold: Optional[float],
) -> None:
    """Add genome-wide significance threshold line.

    Args:
        backend: Plot backend instance.
        ax: Axes object from backend.
        threshold: P-value threshold (e.g., 5e-8). None to skip.
    """
    if threshold is None:
        return
    threshold_line = -np.log10(threshold)
    backend.axhline(
        ax,
        y=threshold_line,
        color=SIGNIFICANCE_LINE_COLOR,
        linestyle="--",
        linewidth=1,
        zorder=1,
    )


def calculate_gene_track_rows(
    genes_df: pd.DataFrame, chrom: int, start: int, end: int
) -> int:
    """Calculate number of gene track rows needed for a region.

    Filters genes to the specified region and computes how many
    vertical rows are needed to display overlapping genes without
    collision.

    Args:
        genes_df: Gene annotations DataFrame with chr, start, end columns.
        chrom: Chromosome number.
        start: Region start position.
        end: Region end position.

    Returns:
        Number of gene rows (minimum 1).
    """
    from .gene_track import assign_gene_positions
    from .utils import normalize_chrom

    chrom_str = normalize_chrom(chrom)
    region_genes = genes_df[
        (genes_df["chr"].astype(str).str.replace("chr", "", regex=False) == chrom_str)
        & (genes_df["end"] >= start)
        & (genes_df["start"] <= end)
    ]
    if not region_genes.empty:
        temp_positions = assign_gene_positions(
            region_genes.sort_values("start"), start, end
        )
        return max(temp_positions) + 1 if temp_positions else 1
    return 1


def calculate_gene_track_height(
    genes_df: pd.DataFrame, chrom: int, start: int, end: int
) -> float:
    """Calculate subplot height-ratio units for a gene track panel.

    Extracted from duplicated logic in LocusZoomPlotter.plot() and
    plot_stacked(). The height grows linearly with the number of stacked
    gene rows so multi-row regions stay legible.

    Args:
        genes_df: Gene annotations DataFrame.
        chrom: Chromosome number.
        start: Region start position.
        end: Region end position.

    Returns:
        Height in the same units as other panel height_ratios.
    """
    n_rows = calculate_gene_track_rows(genes_df, chrom, start, end)
    base_height = 1.0
    per_row_height = 0.5
    return base_height + (n_rows - 1) * per_row_height

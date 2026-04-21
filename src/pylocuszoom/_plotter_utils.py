"""Shared utilities for plotter classes.

Internal module - not part of public API.
"""

from typing import Any, Optional

import numpy as np
import pandas as pd

# Significance thresholds
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8

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


def transform_pvalues(df: pd.DataFrame, p_col: str) -> pd.DataFrame:
    """Validate, filter, and -log10 transform p-values.

    Filters out invalid rows before transformation:
    - NaN p-values are removed.
    - Out-of-range p-values (< 0 or > 1) are removed.
    - Very small valid p-values (< 1e-300) are clipped to avoid -inf.

    Args:
        df: DataFrame with p-value column.
        p_col: Name of p-value column.

    Returns:
        DataFrame with invalid rows removed and neglog10p column added.
    """
    from .logging import logger

    df = df.copy()
    initial_count = len(df)

    # Filter NaN p-values
    nan_mask = df[p_col].isna()
    nan_count = nan_mask.sum()
    if nan_count > 0:
        logger.warning("Found {} NaN p-values, filtering out", nan_count)
        df = df[~nan_mask]

    # Filter out-of-range p-values
    invalid_mask = (df[p_col] < 0) | (df[p_col] > 1)
    invalid_count = invalid_mask.sum()
    if invalid_count > 0:
        logger.warning(
            "Found {} p-values outside [0, 1] range, filtering out",
            invalid_count,
        )
        df = df[~invalid_mask]

    # Log clipped values at debug level
    clipped_count = (df[p_col] < 1e-300).sum()
    if clipped_count > 0:
        logger.debug("Clipping {} p-values below 1e-300 to 1e-300", clipped_count)

    filtered_count = initial_count - len(df)
    if filtered_count > 0:
        logger.debug(
            "P-value filtering removed {} of {} rows",
            filtered_count,
            initial_count,
        )

    df["neglog10p"] = -np.log10(df[p_col].clip(lower=1e-300))
    return df


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

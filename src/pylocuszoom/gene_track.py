"""Gene track data and geometry for regional association plots.

Holds the region filter, the greedy row assignment that stacks overlapping
genes, and the strand-arrow geometry.  The drawing that consumes them lives in
``_regional_panels.GenePanel.draw``.
"""

from typing import List, Optional, Union

import pandas as pd

from .colors import STRAND_ARROW_COLORS
from .utils import normalize_chrom, normalize_chrom_series

# Layout constants
ROW_HEIGHT = 0.35  # Total height per row (reduced for tighter spacing)
GENE_AREA = 0.25  # Bottom portion for gene drawing
EXON_HEIGHT = 0.20  # Exon rectangle height
INTRON_HEIGHT = 0.02  # Thin intron line

# Arrow dimensions (pre-computed for clarity)
ARROW_HEIGHT_RATIO = 0.2625  # EXON_HEIGHT * 0.35 * 0.75 (75% of original height)
ARROW_WIDTH_RATIO = 0.0066  # region_width * 0.006 * 1.1 (10% wider than original)


def assign_gene_positions(genes_df: pd.DataFrame, start: int, end: int) -> List[int]:
    """Assign row indices to genes to minimize overlap.

    Uses a greedy algorithm to stack genes vertically, placing each gene
    in the lowest row where it doesn't overlap with existing genes.

    Args:
        genes_df: Gene annotations DataFrame sorted by start position.
        start: Region start position.
        end: Region end position.

    Returns:
        List of integer row indices (0, 1, 2, ...) for each gene.
    """
    positions = []
    # Track the rightmost end position for each row (including label buffer)
    row_ends: dict[int, int] = {}  # row -> rightmost end position
    region_width = end - start
    label_buffer = region_width * 0.08  # Extra space for labels

    for _, gene in genes_df.iterrows():
        gene_start = max(gene["start"], start)
        gene_end = min(gene["end"], end)

        # Find first available row where gene doesn't overlap
        row = 0
        while row in row_ends and row_ends[row] > gene_start - label_buffer:
            row += 1

        positions.append(row)
        # Update the row's end position (including buffer for next gene check)
        row_ends[row] = gene_end

    return positions


def get_nearest_gene(
    genes_df: pd.DataFrame,
    chrom: Union[int, str],
    pos: int,
    window: int = 50000,
) -> Optional[str]:
    """Get the nearest gene name for a genomic position.

    Searches for genes that overlap or are within the specified window
    of the given position, returning the closest by midpoint distance.

    Args:
        genes_df: Gene annotations DataFrame with chr, start, end, gene_name.
        chrom: Chromosome number or string.
        pos: Position in base pairs.
        window: Window size in bp for searching nearby genes.

    Returns:
        Gene name string or None if no gene found within window.

    Example:
        >>> gene = get_nearest_gene(genes_df, chrom=1, pos=1500000)
        >>> gene
        'BRCA1'
    """
    nearby = filter_genes_by_region(genes_df, chrom, pos - window, pos + window)

    if nearby.empty:
        return None

    # Return the closest gene (by midpoint distance)
    nearby["dist"] = abs((nearby["start"] + nearby["end"]) / 2 - pos)
    return nearby.loc[nearby["dist"].idxmin(), "gene_name"]


def filter_genes_by_region(
    df: pd.DataFrame, chrom: Union[int, str], start: int, end: int
) -> pd.DataFrame:
    """Filter a DataFrame to genes/exons overlapping a genomic region."""
    return df[
        (normalize_chrom_series(df["chr"]) == normalize_chrom(chrom))
        & (df["end"] >= start)
        & (df["start"] <= end)
    ].copy()


def compute_arrow_geometry(
    gene_start: int, gene_end: int, region_width: int, strand: str
) -> tuple[list[float], float, float, str]:
    """Compute arrow tip positions and dimensions for strand arrows.

    Returns:
        Tuple of (arrow_tip_positions, tri_height, tri_width, arrow_color).
    """
    tri_height = EXON_HEIGHT * ARROW_HEIGHT_RATIO
    tri_width = region_width * ARROW_WIDTH_RATIO

    tip_offset = tri_width / 2
    tail_offset = tri_width * 1.5
    gene_center = (gene_start + gene_end) / 2

    if strand == "+":
        arrow_tip_positions = [
            gene_start + tail_offset,
            gene_center + tri_width / 2,
            gene_end - tip_offset,
        ]
        arrow_color = STRAND_ARROW_COLORS["+"]
    else:
        arrow_tip_positions = [
            gene_end - tail_offset,
            gene_center - tri_width / 2,
            gene_start + tip_offset,
        ]
        arrow_color = STRAND_ARROW_COLORS["-"]

    return arrow_tip_positions, tri_height, tri_width, arrow_color

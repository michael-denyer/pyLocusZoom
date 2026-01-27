"""Gene track visualization for regional association plots.

Provides LocusZoom-style gene track plotting with:
- Thin horizontal lines for introns
- Thick rectangles for exons
- Arrows indicating strand direction
- Gene name labels
"""

from typing import Any, List, Optional, Union

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Polygon, Rectangle

from .utils import normalize_chrom

# Strand-specific colors (distinct from LD palette)
STRAND_COLORS: dict[Optional[str], str] = {
    "+": "#DAA520",  # Goldenrod for forward strand
    "-": "#6BB3FF",  # Light blue for reverse strand
    None: "#999999",  # Light grey if no strand info
}

# Layout constants
ROW_HEIGHT = 0.35  # Total height per row (reduced for tighter spacing)
GENE_AREA = 0.25  # Bottom portion for gene drawing
EXON_HEIGHT = 0.20  # Exon rectangle height
INTRON_HEIGHT = 0.02  # Thin intron line


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
    occupied = []  # List of (end_pos, row)
    region_width = end - start

    for _, gene in genes_df.iterrows():
        gene_start = max(gene["start"], start)
        gene_end = min(gene["end"], end)

        # Find first available row with buffer for label spacing
        row = 0
        label_buffer = region_width * 0.08  # Extra space for labels
        for occ_end, occ_row in occupied:
            if occ_row == row and occ_end > gene_start - label_buffer:
                row = occ_row + 1

        positions.append(row)
        occupied.append((gene_end, row))

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
    chrom_str = normalize_chrom(chrom)
    chrom_genes = genes_df[
        genes_df["chr"].astype(str).str.replace("chr", "", regex=False) == chrom_str
    ]

    if chrom_genes.empty:
        return None

    # Find genes that overlap or are within window
    nearby = chrom_genes[
        (chrom_genes["start"] - window <= pos) & (chrom_genes["end"] + window >= pos)
    ]

    if nearby.empty:
        return None

    # Return the closest gene (by midpoint distance)
    nearby = nearby.copy()
    nearby["dist"] = abs((nearby["start"] + nearby["end"]) / 2 - pos)
    return nearby.loc[nearby["dist"].idxmin(), "gene_name"]


def plot_gene_track(
    ax: Axes,
    genes_df: pd.DataFrame,
    chrom: Union[int, str],
    start: int,
    end: int,
    exons_df: Optional[pd.DataFrame] = None,
) -> None:
    """Plot gene annotations as a LocusZoom-style track.

    Creates a gene track with:
    - Thin horizontal lines for introns (gene body)
    - Thick rectangles for exons
    - Arrows indicating strand direction
    - Gene name labels

    Args:
        ax: Matplotlib axes for gene track.
        genes_df: Gene annotations with chr, start, end, gene_name,
            and optionally strand (+/-) column.
        chrom: Chromosome number or string.
        start: Region start position.
        end: Region end position.
        exons_df: Exon annotations with chr, start, end, gene_name
            columns for drawing exon structure. Optional.
    """
    chrom_str = normalize_chrom(chrom)
    region_genes = genes_df[
        (genes_df["chr"].astype(str).str.replace("chr", "", regex=False) == chrom_str)
        & (genes_df["end"] >= start)
        & (genes_df["start"] <= end)
    ].copy()

    ax.set_xlim(start, end)
    ax.set_ylabel("")
    ax.set_yticks([])

    # theme_classic: only bottom spine
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)

    if region_genes.empty:
        ax.set_ylim(0, 1)
        ax.text(
            (start + end) / 2,
            0.5,
            "No genes",
            ha="center",
            va="center",
            fontsize=9,
            color="grey",
            style="italic",
        )
        return

    # Assign vertical positions to avoid overlap
    region_genes = region_genes.sort_values("start")
    positions = assign_gene_positions(region_genes, start, end)

    # Set y-axis limits - small bottom margin for gene body, tight top
    max_row = max(positions) if positions else 0
    bottom_margin = EXON_HEIGHT / 2 + 0.02  # Room for bottom gene
    top_margin = 0.05  # Minimal space above top label
    ax.set_ylim(
        -bottom_margin,
        (max_row + 1) * ROW_HEIGHT - ROW_HEIGHT + GENE_AREA + top_margin,
    )

    # Filter exons for this region if available
    region_exons = None
    if exons_df is not None and not exons_df.empty:
        region_exons = exons_df[
            (
                exons_df["chr"].astype(str).str.replace("chr", "", regex=False)
                == chrom_str
            )
            & (exons_df["end"] >= start)
            & (exons_df["start"] <= end)
        ].copy()

    region_width = end - start

    for idx, (_, gene) in enumerate(region_genes.iterrows()):
        gene_start = max(int(gene["start"]), start)
        gene_end = min(int(gene["end"]), end)
        row = positions[idx]
        gene_name = gene.get("gene_name", "")

        # Get strand-specific color
        strand = gene.get("strand") if "strand" in gene.index else None
        gene_col = STRAND_COLORS.get(strand, STRAND_COLORS[None])

        # Y position: bottom of row + offset for gene area
        y_gene = row * ROW_HEIGHT + 0.05
        y_label = y_gene + EXON_HEIGHT / 2 + 0.01  # Just above gene top

        # Check if we have exon data for this gene
        gene_exons = None
        if region_exons is not None and not region_exons.empty and gene_name:
            gene_exons = region_exons[region_exons["gene_name"] == gene_name].copy()

        if gene_exons is not None and not gene_exons.empty:
            # Draw intron line (thin horizontal line spanning gene)
            ax.add_patch(
                Rectangle(
                    (gene_start, y_gene - INTRON_HEIGHT / 2),
                    gene_end - gene_start,
                    INTRON_HEIGHT,
                    facecolor=gene_col,
                    edgecolor=gene_col,
                    linewidth=0.5,
                    zorder=1,
                )
            )

            # Draw exons (thick rectangles)
            for _, exon in gene_exons.iterrows():
                exon_start = max(int(exon["start"]), start)
                exon_end = min(int(exon["end"]), end)
                ax.add_patch(
                    Rectangle(
                        (exon_start, y_gene - EXON_HEIGHT / 2),
                        exon_end - exon_start,
                        EXON_HEIGHT,
                        facecolor=gene_col,
                        edgecolor=gene_col,
                        linewidth=0.5,
                        zorder=2,
                    )
                )
        else:
            # No exon data - draw full gene body as rectangle (fallback)
            ax.add_patch(
                Rectangle(
                    (gene_start, y_gene - EXON_HEIGHT / 2),
                    gene_end - gene_start,
                    EXON_HEIGHT,
                    facecolor=gene_col,
                    edgecolor=gene_col,
                    linewidth=0.5,
                    zorder=2,
                )
            )

        # Add strand direction triangles (tip, center, tail)
        if "strand" in gene.index:
            strand = gene["strand"]
            arrow_dir = 1 if strand == "+" else -1

            # Triangle dimensions
            tri_height = EXON_HEIGHT * 0.35
            tri_width = region_width * 0.006

            # Arrow positions: front, middle, back (tip positions)
            tip_offset = tri_width / 2  # Tiny offset to keep tip inside gene
            tail_offset = tri_width * 1.5  # Offset for tail arrow from gene start/end
            gene_center = (gene_start + gene_end) / 2
            if arrow_dir == 1:  # Forward strand
                arrow_tip_positions = [
                    gene_start + tail_offset,  # Tail (tip inside gene)
                    gene_center + tri_width / 2,  # Middle (arrow center at gene center)
                    gene_end - tip_offset,  # Tip (near gene end)
                ]
                arrow_color = "#000000"  # Black for forward
            else:  # Reverse strand
                arrow_tip_positions = [
                    gene_end - tail_offset,  # Tail (tip inside gene)
                    gene_center - tri_width / 2,  # Middle (arrow center at gene center)
                    gene_start + tip_offset,  # Tip (near gene start)
                ]
                arrow_color = "#333333"  # Dark grey for reverse

            for tip_x in arrow_tip_positions:
                if arrow_dir == 1:
                    base_x = tip_x - tri_width
                    tri_points = [
                        [tip_x, y_gene],  # Tip pointing right
                        [base_x, y_gene + tri_height],
                        [base_x, y_gene - tri_height],
                    ]
                else:
                    base_x = tip_x + tri_width
                    tri_points = [
                        [tip_x, y_gene],  # Tip pointing left
                        [base_x, y_gene + tri_height],
                        [base_x, y_gene - tri_height],
                    ]

                triangle = Polygon(
                    tri_points,
                    closed=True,
                    facecolor=arrow_color,
                    edgecolor=arrow_color,
                    linewidth=0.5,
                    zorder=5,
                )
                ax.add_patch(triangle)

        # Add gene name label in the gap above gene
        if gene_name:
            label_pos = (gene_start + gene_end) / 2
            ax.text(
                label_pos,
                y_label,
                gene_name,
                ha="center",
                va="bottom",
                fontsize=5.5,
                color="#000000",
                fontweight="medium",
                style="italic",
                zorder=4,
                clip_on=True,
            )


def plot_gene_track_generic(
    ax: Any,
    backend: Any,
    genes_df: pd.DataFrame,
    chrom: Union[int, str],
    start: int,
    end: int,
    exons_df: Optional[pd.DataFrame] = None,
) -> None:
    """Plot gene annotations using a backend-agnostic approach.

    This function works with matplotlib, plotly, and bokeh backends.

    Args:
        ax: Axes object (format depends on backend).
        backend: Backend instance with drawing methods.
        genes_df: Gene annotations with chr, start, end, gene_name,
            and optionally strand (+/-) column.
        chrom: Chromosome number or string.
        start: Region start position.
        end: Region end position.
        exons_df: Exon annotations with chr, start, end, gene_name
            columns for drawing exon structure. Optional.
    """
    chrom_str = normalize_chrom(chrom)
    region_genes = genes_df[
        (genes_df["chr"].astype(str).str.replace("chr", "", regex=False) == chrom_str)
        & (genes_df["end"] >= start)
        & (genes_df["start"] <= end)
    ].copy()

    backend.set_xlim(ax, start, end)
    backend.set_ylabel(ax, "", fontsize=10)
    backend.hide_yaxis(ax)

    if region_genes.empty:
        backend.set_ylim(ax, 0, 1)
        backend.add_text(
            ax,
            (start + end) / 2,
            0.5,
            "No genes",
            fontsize=9,
            ha="center",
            va="center",
            color="grey",
        )
        return

    # Assign vertical positions to avoid overlap
    region_genes = region_genes.sort_values("start")
    positions = assign_gene_positions(region_genes, start, end)

    # Set y-axis limits - small bottom margin for gene body, tight top
    max_row = max(positions) if positions else 0
    bottom_margin = EXON_HEIGHT / 2 + 0.02  # Room for bottom gene
    top_margin = 0.05  # Minimal space above top label
    backend.set_ylim(
        ax,
        -bottom_margin,
        (max_row + 1) * ROW_HEIGHT - ROW_HEIGHT + GENE_AREA + top_margin,
    )

    # Filter exons for this region if available
    region_exons = None
    if exons_df is not None and not exons_df.empty:
        region_exons = exons_df[
            (
                exons_df["chr"].astype(str).str.replace("chr", "", regex=False)
                == chrom_str
            )
            & (exons_df["end"] >= start)
            & (exons_df["start"] <= end)
        ].copy()

    region_width = end - start

    for idx, (_, gene) in enumerate(region_genes.iterrows()):
        gene_start = max(int(gene["start"]), start)
        gene_end = min(int(gene["end"]), end)
        row = positions[idx]
        gene_name = gene.get("gene_name", "")

        # Get strand-specific color
        strand = gene.get("strand") if "strand" in gene.index else None
        gene_col = STRAND_COLORS.get(strand, STRAND_COLORS[None])

        # Y position: bottom of row + offset for gene area
        y_gene = row * ROW_HEIGHT + 0.05
        y_label = y_gene + EXON_HEIGHT / 2 + 0.01  # Just above gene top

        # Check if we have exon data for this gene
        gene_exons = None
        if region_exons is not None and not region_exons.empty and gene_name:
            gene_exons = region_exons[region_exons["gene_name"] == gene_name].copy()

        if gene_exons is not None and not gene_exons.empty:
            # Draw intron line (thin horizontal line spanning gene)
            backend.add_rectangle(
                ax,
                (gene_start, y_gene - INTRON_HEIGHT / 2),
                gene_end - gene_start,
                INTRON_HEIGHT,
                facecolor=gene_col,
                edgecolor=gene_col,
                linewidth=0.5,
                zorder=1,
            )

            # Draw exons (thick rectangles)
            for _, exon in gene_exons.iterrows():
                exon_start = max(int(exon["start"]), start)
                exon_end = min(int(exon["end"]), end)
                backend.add_rectangle(
                    ax,
                    (exon_start, y_gene - EXON_HEIGHT / 2),
                    exon_end - exon_start,
                    EXON_HEIGHT,
                    facecolor=gene_col,
                    edgecolor=gene_col,
                    linewidth=0.5,
                    zorder=2,
                )
        else:
            # No exon data - draw full gene body as rectangle (fallback)
            backend.add_rectangle(
                ax,
                (gene_start, y_gene - EXON_HEIGHT / 2),
                gene_end - gene_start,
                EXON_HEIGHT,
                facecolor=gene_col,
                edgecolor=gene_col,
                linewidth=0.5,
                zorder=2,
            )

        # Add strand direction triangles (tip, center, tail)
        if "strand" in gene.index:
            strand = gene["strand"]
            arrow_dir = 1 if strand == "+" else -1

            # Triangle dimensions
            tri_height = EXON_HEIGHT * 0.35
            tri_width = region_width * 0.006

            # Arrow positions: front, middle, back (tip positions)
            tip_offset = tri_width / 2  # Tiny offset to keep tip inside gene
            tail_offset = tri_width * 1.5  # Offset for tail arrow from gene start/end
            gene_center = (gene_start + gene_end) / 2
            if arrow_dir == 1:  # Forward strand
                arrow_tip_positions = [
                    gene_start + tail_offset,  # Tail (tip inside gene)
                    gene_center + tri_width / 2,  # Middle (arrow center at gene center)
                    gene_end - tip_offset,  # Tip (near gene end)
                ]
                arrow_color = "#000000"  # Black for forward
            else:  # Reverse strand
                arrow_tip_positions = [
                    gene_end - tail_offset,  # Tail (tip inside gene)
                    gene_center - tri_width / 2,  # Middle (arrow center at gene center)
                    gene_start + tip_offset,  # Tip (near gene start)
                ]
                arrow_color = "#333333"  # Dark grey for reverse

            for tip_x in arrow_tip_positions:
                if arrow_dir == 1:
                    base_x = tip_x - tri_width
                    tri_points = [
                        [tip_x, y_gene],  # Tip pointing right
                        [base_x, y_gene + tri_height],
                        [base_x, y_gene - tri_height],
                    ]
                else:
                    base_x = tip_x + tri_width
                    tri_points = [
                        [tip_x, y_gene],  # Tip pointing left
                        [base_x, y_gene + tri_height],
                        [base_x, y_gene - tri_height],
                    ]

                backend.add_polygon(
                    ax,
                    tri_points,
                    facecolor=arrow_color,
                    edgecolor=arrow_color,
                    linewidth=0.5,
                    zorder=5,
                )

        # Add gene name label in the gap above gene
        if gene_name:
            label_pos = (gene_start + gene_end) / 2
            backend.add_text(
                ax,
                label_pos,
                y_label,
                gene_name,
                fontsize=6,
                ha="center",
                va="bottom",
                color="#000000",
            )

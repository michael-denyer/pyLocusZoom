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


def _filter_genes_by_region(
    df: pd.DataFrame, chrom: Union[int, str], start: int, end: int
) -> pd.DataFrame:
    """Filter a DataFrame to genes/exons within a genomic region."""
    chrom_str = normalize_chrom(chrom)
    return df[
        (df["chr"].astype(str).str.replace("chr", "", regex=False) == chrom_str)
        & (df["end"] >= start)
        & (df["start"] <= end)
    ].copy()


def _compute_arrow_geometry(
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
        arrow_color = "#000000"
    else:
        arrow_tip_positions = [
            gene_end - tail_offset,
            gene_center - tri_width / 2,
            gene_start + tip_offset,
        ]
        arrow_color = "#333333"

    return arrow_tip_positions, tri_height, tri_width, arrow_color


def _draw_strand_arrows_matplotlib(
    ax: Axes,
    gene: pd.Series,
    gene_start: int,
    gene_end: int,
    y_gene: float,
    region_width: int,
) -> None:
    """Draw strand direction arrows using matplotlib."""
    strand = gene["strand"]
    arrow_tip_positions, tri_height, tri_width, arrow_color = _compute_arrow_geometry(
        gene_start, gene_end, region_width, strand
    )

    for tip_x in arrow_tip_positions:
        if strand == "+":
            base_x = tip_x - tri_width
        else:
            base_x = tip_x + tri_width

        tri_points = [
            [tip_x, y_gene],
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


def _draw_strand_arrows_generic(
    ax: Any,
    backend: Any,
    gene: pd.Series,
    gene_start: int,
    gene_end: int,
    y_gene: float,
    region_width: int,
) -> None:
    """Draw strand direction arrows using a generic backend."""
    strand = gene["strand"]
    arrow_tip_positions, tri_height, tri_width, arrow_color = _compute_arrow_geometry(
        gene_start, gene_end, region_width, strand
    )

    for tip_x in arrow_tip_positions:
        if strand == "+":
            base_x = tip_x - tri_width
        else:
            base_x = tip_x + tri_width

        tri_points = [
            [tip_x, y_gene],
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
    region_genes = _filter_genes_by_region(genes_df, chrom, start, end)

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
        max_row * ROW_HEIGHT + GENE_AREA + top_margin,
    )

    # Filter exons for this region if available
    region_exons = None
    if exons_df is not None and not exons_df.empty:
        region_exons = _filter_genes_by_region(exons_df, chrom, start, end)

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

        # Add strand direction triangles
        if "strand" in gene.index:
            _draw_strand_arrows_matplotlib(
                ax, gene, gene_start, gene_end, y_gene, region_width
            )

        # Add gene name label in the gap above gene
        if gene_name:
            label_pos = (gene_start + gene_end) / 2
            ax.text(
                label_pos,
                y_label,
                gene_name,
                ha="center",
                va="bottom",
                fontsize=9,
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
    region_genes = _filter_genes_by_region(genes_df, chrom, start, end)

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
        max_row * ROW_HEIGHT + GENE_AREA + top_margin,
    )

    # Filter exons for this region if available
    region_exons = None
    if exons_df is not None and not exons_df.empty:
        region_exons = _filter_genes_by_region(exons_df, chrom, start, end)

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

        # Add strand direction triangles
        if "strand" in gene.index:
            _draw_strand_arrows_generic(
                ax, backend, gene, gene_start, gene_end, y_gene, region_width
            )

        # Add gene name label in the gap above gene
        if gene_name:
            label_pos = (gene_start + gene_end) / 2
            backend.add_text(
                ax,
                label_pos,
                y_label,
                gene_name,
                fontsize=9,
                ha="center",
                va="bottom",
                color="#000000",
            )

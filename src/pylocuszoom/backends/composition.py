"""Backend-neutral figure composition above the rendering seam.

Legend construction lives here as pure functions that produce backend-neutral
``LegendEntry`` specs. Backends implement a single ``add_legend`` primitive that
renders a list of entries natively, so the bin-iteration and swatch policy is
owned once instead of duplicated in every adapter.

The same rule applies to any geometry a backend would otherwise recompute:
``render_recombination_overlay`` drives the secondary-axis primitives, and
``heatmap_highlight_rects`` decides where a SNP highlight is drawn, leaving
each adapter to draw plain rectangles.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, List, Optional, Sequence, Tuple

from ..colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    EQTL_NEGATIVE_BINS,
    EQTL_POSITIVE_BINS,
    LD_BINS,
    LD_NA_COLOR,
    LEAD_SNP_COLOR,
    RECOMB_COLOR,
    EQTLBin,
    LDBin,
    get_credible_set_color,
)

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

    from .base import PlotBackend


@dataclass(frozen=True)
class LegendEntry:
    """One backend-neutral legend entry.

    ``marker`` uses the same vocabulary as ``scatter(marker=...)``: ``"patch"``
    for a filled swatch, or a marker code (``"D"``, ``"^"``, ``"v"``, ``"o"``)
    for a point marker. Backends map these to their native symbols.
    """

    label: str
    color: str
    marker: str = "patch"
    edgecolor: Optional[str] = None


# Mathtext, so matplotlib renders an italic r with a true superscript. The
# interactive backends run legend titles through convert_latex_to_unicode and
# show "r²".
LD_LEGEND_TITLE = r"$r^2$"


def ld_legend_entries(
    ld_bins: Sequence[LDBin] = LD_BINS,
    lead_snp_color: str = LEAD_SNP_COLOR,
) -> List[LegendEntry]:
    """Lead-SNP marker followed by one swatch per LD r^2 bin."""
    entries = [LegendEntry("Lead SNP", lead_snp_color, marker="D")]
    entries.extend(LegendEntry(b.label, b.color, marker="patch") for b in ld_bins)
    return entries


def effect_legend_entries() -> List[LegendEntry]:
    """Colocalization effect-direction categories."""
    return [
        LegendEntry("Same direction", EFFECT_CONGRUENT_COLOR, marker="o"),
        LegendEntry("Opposite direction", EFFECT_INCONGRUENT_COLOR, marker="o"),
        LegendEntry("Missing effect", LD_NA_COLOR, marker="o"),
    ]


def eqtl_legend_entries(
    positive_bins: Sequence[EQTLBin] = EQTL_POSITIVE_BINS,
    negative_bins: Sequence[EQTLBin] = EQTL_NEGATIVE_BINS,
) -> List[LegendEntry]:
    """eQTL effect-size bins: up-triangles positive, down-triangles negative."""
    entries = [LegendEntry(b.label, b.color, marker="^") for b in positive_bins]
    entries.extend(LegendEntry(b.label, b.color, marker="v") for b in negative_bins)
    return entries


def finemapping_legend_entries(
    credible_sets: Sequence[int],
    get_color: Callable[[int], str] = get_credible_set_color,
) -> List[LegendEntry]:
    """One marker per credible set."""
    return [
        LegendEntry(f"CS{cs_id}", get_color(cs_id), marker="o")
        for cs_id in credible_sets
    ]


def render_recombination_overlay(
    backend: "PlotBackend",
    ax: object,
    recomb_df: "pd.DataFrame",
    start: int,
    end: int,
) -> None:
    """Draw a recombination-rate track on a backend's secondary axis.

    Filters ``recomb_df`` to ``[start, end]`` and drives the backend's
    secondary axis. ``create_twin_axis`` returns a handle that is both a
    drawing target for ``fill_between`` and ``line`` and the subject of the
    two ``*_secondary`` axis calls, so the composition is identical across
    backends.
    """
    region = recomb_df[(recomb_df["pos"] >= start) & (recomb_df["pos"] <= end)]
    if region.empty:
        return
    secondary = backend.create_twin_axis(ax)
    backend.fill_between(
        secondary, region["pos"], 0, region["rate"], color=RECOMB_COLOR, alpha=0.15
    )
    backend.line(
        secondary,
        region["pos"],
        region["rate"],
        color=RECOMB_COLOR,
        linewidth=2.5,
        alpha=0.8,
    )
    max_rate = region["rate"].max()
    backend.set_secondary_ylim(secondary, 0, max(max_rate * 1.3, 10))
    backend.set_secondary_ylabel(
        secondary, "Recombination rate (cM/Mb)", color="black", fontsize=9
    )


def mb_tick_positions(
    x_min_bp: float, x_max_bp: float
) -> Tuple[List[float], List[str]]:
    """Tick positions and labels for a genomic axis labelled in megabases.

    Tick spacing widens with the span so a region keeps a readable number of
    ticks, from 0.1 Mb on a sub-megabase region to 5 Mb on a whole chromosome.

    Args:
        x_min_bp: Left edge of the axis in base pairs.
        x_max_bp: Right edge of the axis in base pairs.

    Returns:
        Tuple of (tick positions in base pairs, tick labels in Mb).
    """
    x_min_mb, x_max_mb = x_min_bp / 1e6, x_max_bp / 1e6
    span_mb = x_max_mb - x_min_mb
    for limit, step in ((0.5, 0.1), (2, 0.25), (5, 0.5), (20, 2)):
        if span_mb <= limit:
            tick_step = step
            break
    else:
        tick_step = 5

    first_tick = math.ceil(x_min_mb / tick_step) * tick_step
    tickvals_mb = []
    tick = first_tick
    while tick <= x_max_mb + tick_step / 2:
        tickvals_mb.append(tick)
        tick += tick_step
    return [v * 1e6 for v in tickvals_mb], [f"{v:.2f}" for v in tickvals_mb]


def lower_triangle(matrix: "np.ndarray") -> "np.ndarray":
    """Mask a symmetric matrix's upper triangle, keeping the diagonal.

    An LD matrix is symmetric, so a heatmap shows one half of it. Masking above
    the seam means every backend receives data already in the shape it draws,
    rather than each reimplementing the mask.

    Args:
        matrix: 2D array of values, typically an LD matrix.

    Returns:
        A masked copy whose entries above the diagonal are masked out.
    """
    import numpy as np

    mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
    return np.ma.array(matrix, mask=mask)


def heatmap_highlight_cells(snp_idx: int, n_snps: int) -> List[Tuple[int, int]]:
    """Cells to outline when marking one SNP in a lower-triangular LD heatmap.

    Walks the SNP's row left of the diagonal, then its column below the
    diagonal, so every returned cell lies in the rendered half of the matrix.

    Args:
        snp_idx: Index of the SNP to highlight (0-indexed, must be < n_snps).
        n_snps: Total number of SNPs in the matrix (must be > 0).

    Returns:
        List of ``(x, y)`` matrix index pairs, row cells before column cells.

    Raises:
        ValueError: If snp_idx is out of bounds or n_snps < 1.
    """
    if n_snps < 1 or snp_idx < 0 or snp_idx >= n_snps:
        raise ValueError(f"Invalid snp_idx={snp_idx} for n_snps={n_snps}")
    cells = [(j, snp_idx) for j in range(snp_idx + 1)]
    cells.extend((snp_idx, i) for i in range(snp_idx + 1, n_snps))
    return cells


def cell_edges(coords: Sequence[float]) -> List[Tuple[float, float]]:
    """Left and right boundary of each heatmap cell, at midpoints between centres.

    A heatmap places one cell per coordinate, so a cell reaches halfway to each
    neighbour. The outer cells mirror the gap on their populated side.

    Args:
        coords: Cell centres along one axis, in ascending order.

    Returns:
        One ``(low, high)`` boundary pair per coordinate.
    """
    if len(coords) == 1:
        return [(coords[0] - 0.5, coords[0] + 0.5)]
    mids = [(a + b) / 2 for a, b in zip(coords, coords[1:])]
    first = coords[0] - (mids[0] - coords[0])
    last = coords[-1] + (coords[-1] - mids[-1])
    bounds = [first, *mids, last]
    return list(zip(bounds, bounds[1:]))


def heatmap_highlight_rects(
    snp_idx: int, x_coords: Sequence[float], y_coords: Sequence[float]
) -> List[Tuple[float, float, float, float]]:
    """Outline rectangles marking one SNP in a lower-triangular LD heatmap.

    Args:
        snp_idx: Index of the SNP to highlight (0-indexed, must be < n_snps).
        x_coords: Cell centres along the x-axis, one per SNP.
        y_coords: Cell centres along the y-axis, one per SNP.

    Returns:
        One ``(x, y, width, height)`` rectangle per highlighted cell, in the
        same data coordinates the heatmap was drawn in.

    Raises:
        ValueError: If snp_idx is out of bounds or x_coords is empty.
    """
    cells = heatmap_highlight_cells(snp_idx, len(x_coords))
    x_edges, y_edges = cell_edges(x_coords), cell_edges(y_coords)
    rects = []
    for x, y in cells:
        (x0, x1), (y0, y1) = x_edges[x], y_edges[y]
        rects.append((x0, y0, x1 - x0, y1 - y0))
    return rects

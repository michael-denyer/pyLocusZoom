"""Backend-neutral figure composition above the rendering seam.

Legend construction lives here as pure functions that produce backend-neutral
``LegendEntry`` specs. Backends implement a single ``add_legend`` primitive that
renders a list of entries natively, so the bin-iteration and swatch policy is
owned once instead of duplicated in every adapter.

The same rule applies to any geometry a backend would otherwise recompute:
``render_recombination_overlay`` drives the secondary-axis primitives, and
``heatmap_highlight_cells`` decides which matrix cells a SNP highlight covers,
leaving each adapter to draw them.
"""

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
    EQTLBin,
    LDBin,
    get_credible_set_color,
)

if TYPE_CHECKING:
    from typing import Protocol

    import pandas as pd

    from .base import PlotBackend, SupportsSecondaryAxis

    class _SecondaryAxisBackend(PlotBackend, SupportsSecondaryAxis, Protocol):
        """A backend carrying both the required primitives and a secondary axis.

        ``render_recombination_overlay`` needs the intersection: the
        ``*_secondary`` methods from the optional protocol plus ``hide_spines``
        from the required one. Python has no intersection type, so callers state
        it by narrowing with ``isinstance(backend, SupportsSecondaryAxis)``.
        """


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
    backend: "_SecondaryAxisBackend",
    ax: object,
    recomb_df: "pd.DataFrame",
    start: int,
    end: int,
) -> None:
    """Draw a recombination-rate track on a backend's secondary axis.

    Filters ``recomb_df`` to ``[start, end]`` and drives the backend's
    secondary-axis primitives. ``create_twin_axis`` returns an opaque handle
    that is passed straight back to each ``*_secondary`` primitive, so the
    composition is identical across backends.
    """
    from ..recombination import RECOMB_COLOR

    region = recomb_df[(recomb_df["pos"] >= start) & (recomb_df["pos"] <= end)]
    if region.empty:
        return
    secondary = backend.create_twin_axis(ax)
    backend.fill_between_secondary(
        secondary, region["pos"], 0, region["rate"], color=RECOMB_COLOR, alpha=0.15
    )
    backend.line_secondary(
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
    backend.hide_spines(secondary, ["top"])


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

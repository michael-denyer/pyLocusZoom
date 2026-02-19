"""LD color schemes for regional association plots.

Implements LocusZoom-style coloring based on R² linkage disequilibrium values.
Colors match the locuszoomr R package color scheme.
"""

import math
from typing import List, NamedTuple, Optional


def _is_missing(value: Optional[float]) -> bool:
    """Check if value is None or NaN."""
    return value is None or (isinstance(value, float) and math.isnan(value))


# =============================================================================
# LD Bins (NamedTuple for self-documenting access)
# =============================================================================


class LDBin(NamedTuple):
    """LD bin definition with threshold, label, and color."""

    threshold: float
    label: str
    color: str


LD_BINS: list[LDBin] = [
    LDBin(0.8, "0.8 - 1.0", "#FF0000"),  # red
    LDBin(0.6, "0.6 - 0.8", "#FFA500"),  # orange
    LDBin(0.4, "0.4 - 0.6", "#00CD00"),  # green3
    LDBin(0.2, "0.2 - 0.4", "#00EEEE"),  # cyan2
    LDBin(0.0, "0.0 - 0.2", "#4169E1"),  # royalblue
]

LD_NA_COLOR = "#BEBEBE"  # grey - SNPs lacking LD information
LD_NA_LABEL = "NA"

# Lead SNP color (purple diamond)
LEAD_SNP_COLOR = "#7D26CD"  # purple3

# Fine-mapping/SuSiE credible set colors
# Colors for up to 10 credible sets, matching locuszoomr style
CREDIBLE_SET_COLORS: List[str] = [
    "#FF7F00",  # orange (CS1)
    "#1F78B4",  # blue (CS2)
    "#33A02C",  # green (CS3)
    "#E31A1C",  # red (CS4)
    "#6A3D9A",  # purple (CS5)
    "#B15928",  # brown (CS6)
    "#FB9A99",  # pink (CS7)
    "#A6CEE3",  # light blue (CS8)
    "#B2DF8A",  # light green (CS9)
    "#FDBF6F",  # light orange (CS10)
]

# PIP line color (when not showing credible sets)
PIP_LINE_COLOR = "#FF7F00"  # orange

# =============================================================================
# eQTL Effect Bins (NamedTuple for self-documenting access)
# =============================================================================


class EQTLBin(NamedTuple):
    """eQTL effect bin with range, label, and color."""

    min_val: float
    max_val: float
    label: str
    color: str


# Positive effects (upward triangles)
EQTL_POSITIVE_BINS: list[EQTLBin] = [
    EQTLBin(0.3, 0.4, "0.3 : 0.4", "#8B1A1A"),  # dark red/maroon
    EQTLBin(0.2, 0.3, "0.2 : 0.3", "#FF6600"),  # orange
    EQTLBin(0.1, 0.2, "0.1 : 0.2", "#FFB347"),  # light orange
]
# Negative effects (downward triangles)
EQTL_NEGATIVE_BINS: list[EQTLBin] = [
    EQTLBin(-0.2, -0.1, "-0.2 : -0.1", "#66CDAA"),  # medium aquamarine
    EQTLBin(-0.3, -0.2, "-0.3 : -0.2", "#4682B4"),  # steel blue
    EQTLBin(-0.4, -0.3, "-0.4 : -0.3", "#00008B"),  # dark blue
]


def _find_eqtl_bin(effect: float) -> EQTLBin:
    """Find the eQTL bin for a given effect size.

    Shared lookup used by both get_eqtl_color and get_eqtl_bin to avoid
    duplicating the bin-matching logic.

    Args:
        effect: Effect size (beta coefficient). Must not be None/NaN.

    Returns:
        Matching EQTLBin.
    """
    if effect >= 0:
        for b in EQTL_POSITIVE_BINS:
            if b.min_val <= effect < b.max_val or (
                b.max_val == 0.4 and effect >= b.max_val
            ):
                return b
        return EQTL_POSITIVE_BINS[-1]
    else:
        for b in EQTL_NEGATIVE_BINS:
            if b.min_val < effect <= b.max_val or (
                b.min_val == -0.4 and effect <= b.min_val
            ):
                return b
        return EQTL_NEGATIVE_BINS[-1]


def get_eqtl_color(effect: Optional[float]) -> str:
    """Get color based on eQTL effect size.

    Args:
        effect: Effect size (beta coefficient).

    Returns:
        Hex color code string.
    """
    if _is_missing(effect):
        return LD_NA_COLOR
    return _find_eqtl_bin(effect).color


def get_eqtl_bin(effect: Optional[float]) -> str:
    """Get eQTL effect bin label.

    Args:
        effect: Effect size (beta coefficient).

    Returns:
        Bin label string.
    """
    if _is_missing(effect):
        return LD_NA_LABEL
    return _find_eqtl_bin(effect).label


def get_eqtl_color_palette() -> dict[str, str]:
    """Get color palette for eQTL effect bins.

    Returns:
        Dictionary mapping bin labels to hex colors.
    """
    palette = {}
    for b in EQTL_POSITIVE_BINS:
        palette[b.label] = b.color
    for b in EQTL_NEGATIVE_BINS:
        palette[b.label] = b.color
    return palette


def get_ld_color(r2: Optional[float]) -> str:
    """Get LocusZoom-style color based on LD R² value.

    Uses the locuszoomr R package color scheme:
    - 0.8-1.0: red
    - 0.6-0.8: orange
    - 0.4-0.6: green
    - 0.2-0.4: cyan
    - 0.0-0.2: blue
    - NA: grey

    Args:
        r2: R² value between 0 and 1, or NaN for missing LD.

    Returns:
        Hex color code string.

    Example:
        >>> get_ld_color(0.85)
        '#FF0000'
        >>> get_ld_color(0.5)
        '#00CD00'
        >>> get_ld_color(float('nan'))
        '#BEBEBE'
    """
    if _is_missing(r2):
        return LD_NA_COLOR

    for ld_bin in LD_BINS:
        if r2 >= ld_bin.threshold:
            return ld_bin.color

    return LD_BINS[-1].color


def get_ld_bin(r2: Optional[float]) -> str:
    """Get LD bin label for categorical coloring.

    Args:
        r2: R² value between 0 and 1, or NaN for missing LD.

    Returns:
        Bin label string (e.g., "0.8 - 1.0" or "NA").

    Example:
        >>> get_ld_bin(0.85)
        '0.8 - 1.0'
        >>> get_ld_bin(float('nan'))
        'NA'
    """
    if _is_missing(r2):
        return LD_NA_LABEL

    for ld_bin in LD_BINS:
        if r2 >= ld_bin.threshold:
            return ld_bin.label

    return LD_BINS[-1].label


def get_ld_color_palette() -> dict[str, str]:
    """Get color palette mapping bin labels to colors.

    Returns:
        Dictionary mapping bin labels to hex colors, suitable for
        use with seaborn or matplotlib.

    Example:
        >>> palette = get_ld_color_palette()
        >>> palette["0.8 - 1.0"]
        '#FF0000'
    """
    palette = {ld_bin.label: ld_bin.color for ld_bin in LD_BINS}
    palette[LD_NA_LABEL] = LD_NA_COLOR
    return palette


def get_credible_set_color(cs_id: int) -> str:
    """Get color for a credible set.

    Args:
        cs_id: Credible set ID (1-indexed).

    Returns:
        Hex color code string.

    Example:
        >>> get_credible_set_color(1)
        '#FF7F00'
    """
    if cs_id < 1:
        return LD_NA_COLOR
    # Use modulo to cycle through colors if more than 10 credible sets
    idx = (cs_id - 1) % len(CREDIBLE_SET_COLORS)
    return CREDIBLE_SET_COLORS[idx]


def get_credible_set_color_palette(n_sets: int = 10) -> dict[int, str]:
    """Get color palette for credible sets.

    Args:
        n_sets: Number of credible sets to include.

    Returns:
        Dictionary mapping credible set IDs (1-indexed) to hex colors.

    Example:
        >>> palette = get_credible_set_color_palette(3)
        >>> palette[1]
        '#FF7F00'
    """
    return {
        i + 1: CREDIBLE_SET_COLORS[i % len(CREDIBLE_SET_COLORS)] for i in range(n_sets)
    }


# PheWAS category colors - distinct colors for phenotype categories
PHEWAS_CATEGORY_COLORS: List[str] = [
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#984EA3",  # purple
    "#FF7F00",  # orange
    "#FFFF33",  # yellow
    "#A65628",  # brown
    "#F781BF",  # pink
    "#999999",  # grey
    "#66C2A5",  # teal
    "#FC8D62",  # salmon
    "#8DA0CB",  # periwinkle
]


def get_phewas_category_color(category_idx: int) -> str:
    """Get color for a PheWAS category by index.

    Args:
        category_idx: Zero-indexed category number.

    Returns:
        Hex color code string.
    """
    return PHEWAS_CATEGORY_COLORS[category_idx % len(PHEWAS_CATEGORY_COLORS)]


def get_phewas_category_palette(categories: List[str]) -> dict[str, str]:
    """Get color palette mapping category names to colors.

    Args:
        categories: List of unique category names.

    Returns:
        Dictionary mapping category names to hex colors.
    """
    return {cat: get_phewas_category_color(i) for i, cat in enumerate(categories)}


# =============================================================================
# LD Heatmap Colors
# =============================================================================

# Custom colormap name for LD heatmaps
LD_HEATMAP_CMAP_NAME = "ld_heatmap"

# White-to-red gradient for R² heatmaps (0 = white, 1 = red)
LD_HEATMAP_COLORS: List[str] = ["#FFFFFF", "#FF0000"]

# Color for missing/NaN LD values in heatmaps
LD_HEATMAP_MISSING_COLOR = "#808080"  # grey

# Highlight colors for lead and secondary SNPs in heatmaps
LEAD_SNP_HIGHLIGHT_COLOR = "#FF0000"  # red border/outline for lead SNP
SECONDARY_HIGHLIGHT_COLOR = "#0000FF"  # blue for other highlighted SNPs

# =============================================================================
# Colocalization Effect Direction Colors
# =============================================================================

# Colors for effect direction agreement in colocalization plots
EFFECT_CONGRUENT_COLOR = "#4DAF4A"  # green - same direction
EFFECT_INCONGRUENT_COLOR = "#E41A1C"  # red - opposite direction

"""LD color schemes for regional association plots.

Implements LocusZoom-style coloring based on R² linkage disequilibrium values.
Colors match the locuszoomr R package color scheme.
"""

import math
from typing import List, Optional, Tuple


def _is_missing(value: Optional[float]) -> bool:
    """Check if value is None or NaN."""
    return value is None or (isinstance(value, float) and math.isnan(value))


# LD bin thresholds, labels, and colors
# Format: (threshold, label, color)
LD_BINS: List[Tuple[float, str, str]] = [
    (0.8, "0.8 - 1.0", "#FF0000"),  # red
    (0.6, "0.6 - 0.8", "#FFA500"),  # orange
    (0.4, "0.4 - 0.6", "#00CD00"),  # green3
    (0.2, "0.2 - 0.4", "#00EEEE"),  # cyan2
    (0.0, "0.0 - 0.2", "#4169E1"),  # royalblue
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

# eQTL effect size bins - matches locuszoomr color scheme
# Format: (min_threshold, max_threshold, label, color)
# Positive effects (upward triangles)
EQTL_POSITIVE_BINS: List[Tuple[float, float, str, str]] = [
    (0.3, 0.4, "0.3 : 0.4", "#8B1A1A"),  # dark red/maroon
    (0.2, 0.3, "0.2 : 0.3", "#FF6600"),  # orange
    (0.1, 0.2, "0.1 : 0.2", "#FFB347"),  # light orange
]
# Negative effects (downward triangles)
EQTL_NEGATIVE_BINS: List[Tuple[float, float, str, str]] = [
    (-0.2, -0.1, "-0.2 : -0.1", "#66CDAA"),  # medium aquamarine
    (-0.3, -0.2, "-0.3 : -0.2", "#4682B4"),  # steel blue
    (-0.4, -0.3, "-0.4 : -0.3", "#00008B"),  # dark blue
]


def get_eqtl_color(effect: Optional[float]) -> str:
    """Get color based on eQTL effect size.

    Args:
        effect: Effect size (beta coefficient).

    Returns:
        Hex color code string.
    """
    if _is_missing(effect):
        return LD_NA_COLOR

    if effect >= 0:
        for min_t, max_t, _, color in EQTL_POSITIVE_BINS:
            if min_t <= effect < max_t or (max_t == 0.4 and effect >= max_t):
                return color
        return EQTL_POSITIVE_BINS[-1][3]  # smallest positive bin
    else:
        for min_t, max_t, _, color in EQTL_NEGATIVE_BINS:
            if min_t < effect <= max_t or (min_t == -0.4 and effect <= min_t):
                return color
        return EQTL_NEGATIVE_BINS[-1][3]  # smallest negative bin


def get_eqtl_bin(effect: Optional[float]) -> str:
    """Get eQTL effect bin label.

    Args:
        effect: Effect size (beta coefficient).

    Returns:
        Bin label string.
    """
    if _is_missing(effect):
        return LD_NA_LABEL

    if effect >= 0:
        for min_t, max_t, label, _ in EQTL_POSITIVE_BINS:
            if min_t <= effect < max_t or (max_t == 0.4 and effect >= max_t):
                return label
        return EQTL_POSITIVE_BINS[-1][2]
    else:
        for min_t, max_t, label, _ in EQTL_NEGATIVE_BINS:
            if min_t < effect <= max_t or (min_t == -0.4 and effect <= min_t):
                return label
        return EQTL_NEGATIVE_BINS[-1][2]


def get_eqtl_color_palette() -> dict[str, str]:
    """Get color palette for eQTL effect bins.

    Returns:
        Dictionary mapping bin labels to hex colors.
    """
    palette = {}
    for _, _, label, color in EQTL_POSITIVE_BINS:
        palette[label] = color
    for _, _, label, color in EQTL_NEGATIVE_BINS:
        palette[label] = color
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

    for threshold, _, color in LD_BINS:
        if r2 >= threshold:
            return color

    return LD_BINS[-1][2]


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

    for threshold, label, _ in LD_BINS:
        if r2 >= threshold:
            return label

    return LD_BINS[-1][1]


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
    palette = {label: color for _, label, color in LD_BINS}
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

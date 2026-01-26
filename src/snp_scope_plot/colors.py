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

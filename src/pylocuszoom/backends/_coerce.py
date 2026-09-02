"""Matplotlib-vocabulary coercions shared by the interactive backends.

``PlotBackend`` speaks matplotlib's vocabulary: figure sizes in inches, scatter
sizes as marker area, colours as either one value or one per point. Plotly and
bokeh each need the same translations out of it, so they live here as pure
functions instead of once per adapter.
"""

from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# Web backends render at a nominal 100 dpi, so an inch of figsize is 100 pixels.
_PIXELS_PER_INCH = 100


def pixels(figsize: Tuple[float, float]) -> Tuple[int, int]:
    """Convert a matplotlib figure size in inches to pixels.

    Args:
        figsize: Figure size as (width, height) in inches.

    Returns:
        Tuple of (width, height) in pixels.
    """
    return (
        int(figsize[0] * _PIXELS_PER_INCH),
        int(figsize[1] * _PIXELS_PER_INCH),
    )


def normalize_ratios(ratios: Optional[Sequence[float]]) -> Optional[List[float]]:
    """Scale panel ratios to sum to one, leaving ``None`` alone.

    Args:
        ratios: Relative panel sizes, or None when the caller wants the default
            even split.

    Returns:
        The normalized ratios, or None if none were given.
    """
    if ratios is None:
        return None
    total = sum(ratios)
    return [r / total for r in ratios]


def split_pixels(total: int, ratios: Optional[Sequence[float]], n: int) -> List[int]:
    """Divide a pixel extent into ``n`` parts, by ratio or evenly.

    Args:
        total: The extent to divide, in pixels.
        ratios: Relative sizes, or None for an even split.
        n: Number of parts when there are no ratios.

    Returns:
        One pixel size per part.
    """
    if ratios is None:
        return [total // n] * n
    denominator = sum(ratios)
    return [int(total * r / denominator) for r in ratios]


def marker_diameter(
    sizes: Union[float, Sequence[float], pd.Series],
) -> Union[float, List[float]]:
    """Convert matplotlib marker areas to the diameters web backends expect.

    Matplotlib sizes a marker by area; plotly and bokeh size it by diameter.
    A floor of 6 keeps the smallest markers clickable. A scalar stays scalar,
    because both backends serialize one value more compactly than a repeated
    list, and an exported document carries that difference.

    Args:
        sizes: One area for every point, or a single area for all of them.

    Returns:
        One diameter, or one per point.
    """
    if isinstance(sizes, (int, float)):
        return max(6, sizes**0.5)
    return [max(6, s**0.5) for s in sizes]


def marker_colors(
    colors: Union[str, Sequence[str], pd.Series],
) -> Union[str, List[str]]:
    """Pass a single colour or a per-point sequence through as plain values.

    A scalar stays scalar for the same serialization reason as
    ``marker_diameter``.

    Args:
        colors: One colour for every point, or a single colour for all of them.

    Returns:
        One colour, or one per point.
    """
    if isinstance(colors, str):
        return colors
    if isinstance(colors, (pd.Series, np.ndarray)):
        return list(colors)
    return colors


def per_point(value: Union[str, float, Sequence[Any], pd.Series], n: int) -> List[Any]:
    """Expand a scalar to one value per point, or list a per-point sequence.

    For backends that address every point through a columnar data source, where
    a single value is not accepted in place of a column.

    Args:
        value: One value for every point, or a single value for all of them.
        n: Number of points.

    Returns:
        A list of exactly ``n`` values.
    """
    if isinstance(value, (str, int, float)):
        return [value] * n
    return list(value)


def broadcast(
    value: Union[float, pd.Series, Sequence[Any]], n: int
) -> Union[List[Any], np.ndarray]:
    """Repeat a scalar fill bound across ``n`` points, or pass a sequence through.

    A pandas Series is handed over as its ndarray rather than a list. Bokeh packs
    an ndarray into the document as base64 and a list as plain JSON, so
    materializing it would inflate every exported HTML carrying a filled band.

    Args:
        value: One bound for every point, or a single bound for all of them.
        n: Number of points, used to broadcast a scalar.

    Returns:
        The per-point bounds, as an ndarray when the input was a Series.
    """
    if isinstance(value, (int, float)):
        return [value] * n
    return value.values if isinstance(value, pd.Series) else list(value)

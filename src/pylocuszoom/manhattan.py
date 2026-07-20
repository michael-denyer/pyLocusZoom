"""Manhattan plot data preparation and chromosome ordering."""

from typing import Literal

import colorcet as cc
import numpy as np
import pandas as pd

from .logging import logger

# Species aliases
SPECIES_ALIASES: dict[str, str] = {
    "dog": "canine",
    "cat": "feline",
}

# Chromosome orders for supported species
CHROMOSOME_ORDERS: dict[str, list[str]] = {
    "canine": [str(i) for i in range(1, 39)] + ["X", "Y", "MT"],
    "feline": [
        "A1",
        "A2",
        "A3",
        "B1",
        "B2",
        "B3",
        "B4",
        "C1",
        "C2",
        "D1",
        "D2",
        "D3",
        "D4",
        "E1",
        "E2",
        "E3",
        "X",
        "Y",
        "MT",
    ],
    "human": [str(i) for i in range(1, 23)] + ["X", "Y", "MT"],
}


def get_chromosome_order(
    species: Literal["canine", "feline", "human", "dog", "cat"] | None = None,
    custom_order: list[str] | None = None,
) -> list[str]:
    """Get chromosome order for a species.

    Args:
        species: Species name for built-in order. Supports aliases:
            'dog' -> 'canine', 'cat' -> 'feline'.
        custom_order: Custom chromosome order (overrides species).

    Returns:
        List of chromosome names in display order.

    Raises:
        ValueError: If neither species nor custom_order provided,
            or if species is unknown.
    """
    if custom_order is not None:
        return custom_order
    if species is not None:
        # Resolve aliases
        resolved_species = SPECIES_ALIASES.get(species, species)
        if resolved_species not in CHROMOSOME_ORDERS:
            raise ValueError(
                f"Unknown species '{species}'. "
                f"Use one of {list(CHROMOSOME_ORDERS.keys())} "
                f"(or aliases: {list(SPECIES_ALIASES.keys())}) "
                f"or provide custom_order."
            )
        return CHROMOSOME_ORDERS[resolved_species]
    raise ValueError("Must provide either species or custom_order")


def get_chromosome_colors(n_chromosomes: int) -> list[str]:
    """Get perceptually distinct colors for chromosomes.

    Uses colorcet glasbey_dark palette for good visual
    separation with saturated colors.

    Args:
        n_chromosomes: Number of chromosomes to color.

    Returns:
        List of hex color strings.
    """
    palette = cc.b_glasbey_bw_minc_20_maxl_70
    return [palette[i % len(palette)] for i in range(n_chromosomes)]


def _filter_invalid_pvalues(df: pd.DataFrame, p_col: str) -> pd.DataFrame:
    """Filter rows with invalid p-values and compute -log10(p).

    Drops rows where p is NaN, negative, or > 1, then clips remaining
    values at 1e-300 before -log10 transformation. Always returns a
    copy — the input DataFrame is never mutated.

    Args:
        df: DataFrame with p-value column.
        p_col: Name of p-value column.

    Returns:
        DataFrame with invalid rows removed and _neg_log_p column added.

    Raises:
        ValueError: If all rows have invalid p-values.
    """
    valid_mask = df[p_col].between(0, 1, inclusive="both")
    n_invalid = (~valid_mask).sum()
    if n_invalid > 0:
        logger.warning(f"Dropped {n_invalid} invalid p-values (NaN or outside [0, 1])")

    result = df[valid_mask].copy()

    if len(result) == 0:
        raise ValueError(
            f"All rows have invalid p-values in column '{p_col}' "
            f"(NaN, negative, or > 1). Cannot create plot."
        )

    result["_neg_log_p"] = -np.log10(result[p_col].clip(lower=1e-300))
    return result


def prepare_manhattan_data(
    df: pd.DataFrame,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    p_col: str = "p",
    species: Literal["canine", "feline", "human", "dog", "cat"] | None = None,
    custom_order: list[str] | None = None,
) -> pd.DataFrame:
    """Prepare DataFrame for Manhattan plot rendering.

    Computes cumulative positions for x-axis and assigns chromosome colors.

    Args:
        df: GWAS results DataFrame.
        chrom_col: Column name for chromosome.
        pos_col: Column name for position.
        p_col: Column name for p-value.
        species: Species for chromosome ordering.
        custom_order: Custom chromosome order.

    Returns:
        DataFrame with additional columns:
        - _chrom_str: String-normalized chromosome name
        - _chrom_idx: Integer index for chromosome
        - _cumulative_pos: X-axis position
        - _neg_log_p: -log10(p-value)
        - _color: Hex color for chromosome
    """
    # Validate required columns
    for col, name in [(chrom_col, "chrom"), (pos_col, "pos"), (p_col, "p")]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame (for {name})")

    # Get chromosome order
    chrom_order = get_chromosome_order(species, custom_order)

    # Filter invalid p-values early (returns a copy)
    result = _filter_invalid_pvalues(df, p_col)

    # Normalize chromosome names (handle int vs str)
    result["_chrom_str"] = result[chrom_col].astype(str)

    # Map chromosomes to order index (len(chrom_order) for unknown)
    chrom_to_idx = {chrom: i for i, chrom in enumerate(chrom_order)}
    result["_chrom_idx"] = result["_chrom_str"].map(
        lambda x: chrom_to_idx.get(x, len(chrom_order))
    )

    # Sort by chromosome index then position
    result = result.sort_values(["_chrom_idx", pos_col])

    # Calculate cumulative positions. Single groupby pass replaces repeated
    # boolean masking per chromosome, which was O(N·C) on GWAS-scale data.
    max_by_chrom = result.groupby("_chrom_str", sort=False)[pos_col].max()
    chrom_offsets = {}
    cumulative = 0
    for chrom in chrom_order:
        if chrom in max_by_chrom.index:
            chrom_offsets[chrom] = cumulative
            cumulative += int(max_by_chrom[chrom]) + 1_000_000  # 1Mb gap

    # Handle chromosomes not in order
    unknown_chroms = sorted(set(max_by_chrom.index) - set(chrom_order))
    for chrom in unknown_chroms:
        chrom_offsets[chrom] = cumulative
        cumulative += int(max_by_chrom[chrom]) + 1_000_000

    # Calculate cumulative position (vectorized — avoids apply(axis=1))
    result["_cumulative_pos"] = (
        result["_chrom_str"].map(chrom_offsets).fillna(0) + result[pos_col]
    )

    # Assign colors
    all_chroms = chrom_order + sorted(unknown_chroms)
    colors = get_chromosome_colors(len(all_chroms))
    chrom_to_color = {chrom: colors[i] for i, chrom in enumerate(all_chroms)}
    result["_color"] = result["_chrom_str"].map(chrom_to_color)

    # Calculate chromosome centers for x-axis labels. Single groupby.mean()
    # replaces a per-chromosome boolean mask (O(N·C) -> O(N)).
    centers_by_chrom = result.groupby("_chrom_str", sort=False)[
        "_cumulative_pos"
    ].mean()
    chrom_centers = {
        chrom: centers_by_chrom[chrom]
        for chrom in all_chroms
        if chrom in centers_by_chrom.index
    }

    result.attrs["chrom_centers"] = chrom_centers
    result.attrs["chrom_order"] = all_chroms

    return result


def prepare_categorical_data(
    df: pd.DataFrame,
    category_col: str,
    p_col: str = "p",
    category_order: list[str] | None = None,
) -> pd.DataFrame:
    """Prepare DataFrame for categorical Manhattan plot (PheWAS-style).

    Args:
        df: Results DataFrame with categories and p-values.
        category_col: Column name for category.
        p_col: Column name for p-value.
        category_order: Custom category order.

    Returns:
        DataFrame with additional columns for plotting.
    """
    # Validate required columns
    if category_col not in df.columns:
        raise ValueError(f"Column '{category_col}' not found in DataFrame")
    if p_col not in df.columns:
        raise ValueError(f"Column '{p_col}' not found in DataFrame")

    # Filter invalid p-values early, before derived calculations (returns a copy)
    result = _filter_invalid_pvalues(df, p_col)

    # Get category order
    if category_order is None:
        # Get unique values, drop NaN, convert to strings for consistent sorting
        unique_vals = result[category_col].dropna().unique()
        # Convert all to strings and sort to handle mixed types safely
        category_order = sorted([str(v) for v in unique_vals])

    # Convert category column to string for consistent handling
    result["_cat_str"] = result[category_col].astype(str)

    # Map categories to index (use string values for lookup)
    cat_to_idx = {cat: i for i, cat in enumerate(category_order)}
    result["_cat_idx"] = result["_cat_str"].map(
        lambda x: cat_to_idx.get(x, len(category_order))
    )

    # Use category index as x position (with jitter for multiple points per category)
    rng = np.random.default_rng(42)  # Local RNG for reproducible jitter
    result["_x_pos"] = result["_cat_idx"] + rng.uniform(-0.3, 0.3, size=len(result))

    # Assign colors (use string values for lookup)
    colors = get_chromosome_colors(len(category_order))
    cat_to_color = {cat: colors[i] for i, cat in enumerate(category_order)}
    result["_color"] = result["_cat_str"].map(cat_to_color)

    result.attrs["category_order"] = category_order
    result.attrs["category_centers"] = {cat: i for i, cat in enumerate(category_order)}

    return result

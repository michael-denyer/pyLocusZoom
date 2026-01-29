"""Manhattan plot data preparation and chromosome ordering."""

from typing import Literal

import colorcet as cc
import numpy as np
import pandas as pd

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

    Uses colorcet b_glasbey_bw_minc_20_minl_30 palette for good visual
    separation with mid-range luminance (visible on both light and dark
    backgrounds).

    Args:
        n_chromosomes: Number of chromosomes to color.

    Returns:
        List of hex color strings.
    """
    palette = cc.b_glasbey_bw_minc_20_minl_30
    return [palette[i % len(palette)] for i in range(n_chromosomes)]


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

    # Create working copy
    result = df.copy()

    # Normalize chromosome names (handle int vs str)
    result["_chrom_str"] = result[chrom_col].astype(str)

    # Map chromosomes to order index (-1 for unknown)
    chrom_to_idx = {chrom: i for i, chrom in enumerate(chrom_order)}
    result["_chrom_idx"] = result["_chrom_str"].map(
        lambda x: chrom_to_idx.get(x, len(chrom_order))
    )

    # Sort by chromosome index then position
    result = result.sort_values(["_chrom_idx", pos_col])

    # Calculate cumulative positions
    # First get max position per chromosome
    chrom_offsets = {}
    cumulative = 0
    for chrom in chrom_order:
        chrom_data = result[result["_chrom_str"] == chrom]
        if len(chrom_data) > 0:
            chrom_offsets[chrom] = cumulative
            cumulative += chrom_data[pos_col].max() + 1_000_000  # 1Mb gap

    # Handle chromosomes not in order
    unknown_chroms = set(result["_chrom_str"]) - set(chrom_order)
    for chrom in sorted(unknown_chroms):
        chrom_data = result[result["_chrom_str"] == chrom]
        if len(chrom_data) > 0:
            chrom_offsets[chrom] = cumulative
            cumulative += chrom_data[pos_col].max() + 1_000_000

    # Calculate cumulative position
    result["_cumulative_pos"] = result.apply(
        lambda row: chrom_offsets.get(row["_chrom_str"], 0) + row[pos_col], axis=1
    )

    # Calculate -log10(p)
    result["_neg_log_p"] = -np.log10(result[p_col].clip(lower=1e-300))

    # Assign colors
    all_chroms = chrom_order + sorted(unknown_chroms)
    colors = get_chromosome_colors(len(all_chroms))
    chrom_to_color = {chrom: colors[i] for i, chrom in enumerate(all_chroms)}
    result["_color"] = result["_chrom_str"].map(chrom_to_color)

    # Calculate chromosome centers for x-axis labels
    chrom_centers = {}
    for chrom in all_chroms:
        chrom_data = result[result["_chrom_str"] == chrom]
        if len(chrom_data) > 0:
            chrom_centers[chrom] = chrom_data["_cumulative_pos"].mean()

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

    result = df.copy()

    # Get category order
    if category_order is None:
        category_order = sorted(result[category_col].unique())

    # Map categories to index
    cat_to_idx = {cat: i for i, cat in enumerate(category_order)}
    result["_cat_idx"] = result[category_col].map(
        lambda x: cat_to_idx.get(x, len(category_order))
    )

    # Use category index as x position (with jitter for multiple points per category)
    np.random.seed(42)  # Reproducible jitter
    result["_x_pos"] = result["_cat_idx"] + np.random.uniform(
        -0.3, 0.3, size=len(result)
    )

    # Calculate -log10(p)
    result["_neg_log_p"] = -np.log10(result[p_col].clip(lower=1e-300))

    # Assign colors
    colors = get_chromosome_colors(len(category_order))
    cat_to_color = {cat: colors[i] for i, cat in enumerate(category_order)}
    result["_color"] = result[category_col].map(cat_to_color)

    result.attrs["category_order"] = category_order
    result.attrs["category_centers"] = {cat: i for i, cat in enumerate(category_order)}

    return result

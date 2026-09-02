"""Manhattan plot data preparation and chromosome ordering."""

from dataclasses import dataclass
from typing import Literal, Mapping, Sequence, Tuple, Union

import colorcet as cc
import numpy as np
import pandas as pd

from ._data import prepare_pvalue_data

CHROMOSOME_GAP = 1_000_000

ALL_PVALUES_INVALID = (
    "All rows have invalid p-values in column '{p_col}' "
    "(NaN, negative, or > 1). Cannot create plot."
)

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


@dataclass(frozen=True)
class GenomeLayout:
    """Where each chromosome sits on a shared Manhattan x axis.

    Every frame drawn against one layout puts a given ``(chrom, pos)`` at the
    same x, which is what lets a Miami or stacked figure share x limits and one
    set of ticks.

    Attributes:
        order: Chromosomes in display order, including any the frames do not
            carry. A chromosome's colour is its position in this list.
        offsets: X coordinate of each chromosome's first base. Only the
            chromosomes the frames carry appear.
        colors: Hex colour per chromosome in ``order``.
        centers: Mean x of the points on each chromosome, its tick position.
        x_limits: Padded x span covering every point in every frame.
        total_length: X coordinate one gap past the last chromosome's end.
    """

    order: Tuple[str, ...]
    offsets: Mapping[str, int]
    colors: Mapping[str, str]
    centers: Mapping[str, float]
    x_limits: Tuple[float, float]
    total_length: int

    @property
    def tick_labels(self) -> list[str]:
        """Return the label of every chromosome that carries data."""
        return [chrom for chrom in self.order if chrom in self.centers]

    @property
    def tick_positions(self) -> list[float]:
        """Return the tick position of every chromosome that carries data."""
        return [self.centers[chrom] for chrom in self.tick_labels]

    @classmethod
    def from_frames(
        cls,
        frames: Sequence[pd.DataFrame],
        *,
        chrom_col: str,
        pos_col: str,
        order: Sequence[str],
    ) -> "GenomeLayout":
        """Lay out the genome axis once for every frame that shares it.

        Args:
            frames: Frames already filtered to plottable rows.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            order: Display order of the known chromosomes. Chromosomes the
                frames carry but this order omits are appended sorted.

        Returns:
            The shared layout.
        """
        pooled = pd.DataFrame(
            {
                "_chrom_str": pd.concat(
                    [frame[chrom_col].astype(str) for frame in frames],
                    ignore_index=True,
                ),
                "_pos": pd.concat(
                    [frame[pos_col] for frame in frames], ignore_index=True
                ),
            }
        )
        max_by_chrom = pooled.groupby("_chrom_str", sort=False)["_pos"].max()

        offsets: dict[str, int] = {}
        cumulative = 0
        unknown = sorted(set(max_by_chrom.index) - set(order))
        for chrom in list(order) + unknown:
            if chrom in max_by_chrom.index:
                offsets[chrom] = cumulative
                cumulative += int(max_by_chrom[chrom]) + CHROMOSOME_GAP

        full_order = tuple(order) + tuple(unknown)
        chrom_to_idx = {chrom: i for i, chrom in enumerate(full_order)}
        pooled["_chrom_idx"] = pooled["_chrom_str"].map(
            lambda x: chrom_to_idx.get(x, len(full_order))
        )
        pooled = pooled.sort_values(["_chrom_idx", "_pos"])
        pooled["_cumulative_pos"] = pooled["_chrom_str"].map(offsets) + pooled["_pos"]

        x_min = pooled["_cumulative_pos"].min()
        x_max = pooled["_cumulative_pos"].max()
        padding = (x_max - x_min) * 0.01
        return cls(
            order=full_order,
            offsets=offsets,
            colors=dict(zip(full_order, get_chromosome_colors(len(full_order)))),
            centers=pooled.groupby("_chrom_str", sort=False)["_cumulative_pos"]
            .mean()
            .to_dict(),
            x_limits=(x_min - padding, x_max + padding),
            total_length=cumulative,
        )


@dataclass(frozen=True)
class CategoryLayout:
    """Where each category sits on a categorical (PheWAS-style) x axis.

    Attributes:
        order: Categories in display order.
        colors: Hex colour per category.
    """

    order: Tuple[str, ...]
    colors: Mapping[str, str]

    @property
    def x_limits(self) -> Tuple[float, float]:
        """Return the x span, half a slot either side of the categories."""
        return (-0.5, len(self.order) - 0.5)

    @property
    def tick_labels(self) -> list[str]:
        """Return one label per category."""
        return list(self.order)

    @property
    def tick_positions(self) -> list[float]:
        """Return one tick per category, at its slot index."""
        return [float(index) for index in range(len(self.order))]


PanelLayout = Union[GenomeLayout, CategoryLayout]


def prepare_manhattan_frames(
    dfs: Sequence[pd.DataFrame],
    *,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    p_col: str = "p",
    species: Literal["canine", "feline", "human", "dog", "cat"] | None = None,
    custom_order: list[str] | None = None,
) -> list[pd.DataFrame]:
    """Prepare several GWAS frames against one shared genome layout.

    Every returned frame carries the same layout in ``attrs["layout"]``, so a
    given ``(chrom, pos)`` lands at the same x in all of them.

    Args:
        dfs: GWAS results DataFrames, in panel order.
        chrom_col: Column name for chromosome.
        pos_col: Column name for position.
        p_col: Column name for p-value.
        species: Species for chromosome ordering.
        custom_order: Custom chromosome order.

    Returns:
        One prepared frame per input, in the same order.

    Raises:
        ValueError: If a required column is missing, or if neither species nor
            custom_order names a chromosome order.
    """
    for df in dfs:
        for col, name in [(chrom_col, "chrom"), (pos_col, "pos"), (p_col, "p")]:
            if col not in df.columns:
                raise ValueError(f"Column '{col}' not found in DataFrame (for {name})")

    order = get_chromosome_order(species, custom_order)
    filtered = [
        prepare_pvalue_data(df, p_col, on_empty=ALL_PVALUES_INVALID.format(p_col=p_col))
        for df in dfs
    ]
    layout = GenomeLayout.from_frames(
        filtered, chrom_col=chrom_col, pos_col=pos_col, order=order
    )
    return [
        _apply_genome_layout(frame, chrom_col, pos_col, layout) for frame in filtered
    ]


def _apply_genome_layout(
    result: pd.DataFrame,
    chrom_col: str,
    pos_col: str,
    layout: GenomeLayout,
) -> pd.DataFrame:
    """Place one filtered frame on the shared genome axis."""
    result["_chrom_str"] = result[chrom_col].astype(str)
    chrom_to_idx = {chrom: i for i, chrom in enumerate(layout.order)}
    result["_chrom_idx"] = result["_chrom_str"].map(
        lambda x: chrom_to_idx.get(x, len(layout.order))
    )
    result = result.sort_values(["_chrom_idx", pos_col])
    result["_cumulative_pos"] = (
        result["_chrom_str"].map(layout.offsets) + result[pos_col]
    )
    result["_color"] = result["_chrom_str"].map(layout.colors)
    result.attrs["layout"] = layout
    return result


def prepare_manhattan_data(
    df: pd.DataFrame,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    p_col: str = "p",
    species: Literal["canine", "feline", "human", "dog", "cat"] | None = None,
    custom_order: list[str] | None = None,
) -> pd.DataFrame:
    """Prepare one DataFrame for Manhattan plot rendering.

    The single-frame case of :func:`prepare_manhattan_frames`.

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
        - neglog10p: -log10(p-value)
        - _color: Hex color for chromosome

        The shared :class:`GenomeLayout` is in ``attrs["layout"]``.
    """
    return prepare_manhattan_frames(
        [df],
        chrom_col=chrom_col,
        pos_col=pos_col,
        p_col=p_col,
        species=species,
        custom_order=custom_order,
    )[0]


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

    result = prepare_pvalue_data(
        df, p_col, on_empty=ALL_PVALUES_INVALID.format(p_col=p_col)
    )

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

    layout = CategoryLayout(
        order=tuple(category_order),
        colors=dict(zip(category_order, get_chromosome_colors(len(category_order)))),
    )
    result["_color"] = result["_cat_str"].map(layout.colors)
    result.attrs["layout"] = layout

    return result

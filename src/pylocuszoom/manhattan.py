"""Manhattan plot data preparation and chromosome ordering."""

from dataclasses import dataclass
from typing import Mapping, Sequence, Tuple, Union

import colorcet as cc
import numpy as np
import pandas as pd

from ._data import prepare_pvalue_data
from .config import GenomeWideConfig, resolve_deprecated_columns
from .exceptions import ValidationError
from .schemas import Canonical, validate_gwas_df
from .species import Species, resolve_species

CHROMOSOME_GAP = 1_000_000

ALL_PVALUES_INVALID = (
    "All rows have invalid p-values in column '{p_col}' "
    "(NaN, negative, or > 1). Cannot create plot."
)


def get_chromosome_order(
    species: str | Species | None = None,
    custom_order: list[str] | None = None,
) -> list[str]:
    """Get chromosome order for a species.

    Args:
        species: Species name, alias or record carrying a built-in order.
        custom_order: Custom chromosome order (overrides species).

    Returns:
        List of chromosome names in display order.

    Raises:
        ValidationError: If the species is unknown, or is known but has no
            built-in chromosome order.
        ValueError: If neither species nor custom_order provided.
    """
    if custom_order is not None:
        return custom_order
    record = resolve_species(species)
    if record is None:
        raise ValueError("Must provide either species or custom_order")
    if not record.chromosomes:
        raise ValidationError(
            f"No built-in chromosome order for species {record.key!r}; "
            f"pass custom_order."
        )
    return list(record.chromosomes)


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


@dataclass(frozen=True)
class PreparedManhattan:
    """One frame laid out for a Manhattan-style panel, with its layout.

    ``frame`` carries ``neglog10p``, ``_color`` and the x column the layout
    places it on. ``layout`` is the one every frame prepared in the same call
    shares, so a given point lands at the same x in all of them.
    """

    frame: pd.DataFrame
    layout: PanelLayout


def prepare_genomewide_frames(
    dfs: Sequence[pd.DataFrame],
    config: GenomeWideConfig,
    *,
    species: Species | None,
    rs_col: str | None = None,
) -> list[PreparedManhattan]:
    """Validate each frame against ``config`` and lay them out on one genome.

    The boundary for the genome-wide families: every frame is checked for
    the chromosome, position and p-value columns the config names (and
    ``rs_col`` when given) before any of them is laid out. A frame still
    carrying the pre-4.0 column names is accepted with a deprecation warning.

    Args:
        dfs: GWAS results DataFrames, in panel order.
        config: Column names and chromosome order.
        species: Species whose chromosome order lays the axis out, unless
            ``config.custom_chrom_order`` overrides it.
        rs_col: SNP id column to require as well, or None.

    Raises:
        ValidationError: If a frame is empty or lacks a named column.
    """
    config = resolve_deprecated_columns(dfs[0], config) if dfs else config
    for df in dfs:
        validate_gwas_df(
            df,
            pos_col=config.pos_col,
            p_col=config.p_col,
            rs_col=rs_col,
            chrom_col=config.chrom_col,
        )
    return prepare_manhattan_frames(
        dfs,
        chrom_col=config.chrom_col,
        pos_col=config.pos_col,
        p_col=config.p_col,
        species=species,
        custom_order=config.custom_chrom_order,
    )


def prepare_manhattan_frames(
    dfs: Sequence[pd.DataFrame],
    *,
    chrom_col: str = Canonical.CHROM,
    pos_col: str = Canonical.POS,
    p_col: str = Canonical.P,
    species: str | Species | None = None,
    custom_order: list[str] | None = None,
) -> list[PreparedManhattan]:
    """Prepare several GWAS frames against one shared genome layout.

    Every returned value carries the same :class:`GenomeLayout`, so a given
    ``(chrom, pos)`` lands at the same x in all of them. Pass a one-element
    list for a single panel.

    Args:
        dfs: GWAS results DataFrames, in panel order.
        chrom_col: Column name for chromosome.
        pos_col: Column name for position.
        p_col: Column name for p-value.
        species: Species for chromosome ordering.
        custom_order: Custom chromosome order.

    Returns:
        One prepared value per input, in the same order, each carrying the
        columns ``_chrom_str``, ``_chrom_idx``, ``_cumulative_pos``,
        ``neglog10p`` and ``_color``.

    Raises:
        ValueError: If a required column is missing, or if neither species nor
            custom_order names a chromosome order.
    """
    for df in dfs:
        for col, name in [
            (chrom_col, "chromosome"),
            (pos_col, "position"),
            (p_col, "p-value"),
        ]:
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
        PreparedManhattan(
            _apply_genome_layout(frame, chrom_col, pos_col, layout), layout
        )
        for frame in filtered
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
    return result


def prepare_categorical_data(
    df: pd.DataFrame,
    category_col: str,
    p_col: str = Canonical.P,
    category_order: list[str] | None = None,
) -> PreparedManhattan:
    """Prepare DataFrame for categorical Manhattan plot (PheWAS-style).

    Args:
        df: Results DataFrame with categories and p-values.
        category_col: Column name for category.
        p_col: Column name for p-value.
        category_order: Custom category order.

    Returns:
        The frame with ``_cat_str``, ``_cat_idx``, ``_x_pos``, ``neglog10p``
        and ``_color`` columns, and its :class:`CategoryLayout`.
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
    return PreparedManhattan(result, layout)

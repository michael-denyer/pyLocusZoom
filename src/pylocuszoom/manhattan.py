"""Manhattan plot data preparation and chromosome ordering."""

from dataclasses import dataclass
from typing import Mapping, Sequence, Tuple, Union

import colorcet as cc
import numpy as np
import pandas as pd

from ._data import prepare_pvalue_data
from ._plotter_utils import CHROMOSOME_GAP
from .config import GenomeWideConfig, GenomeWideStyle
from .exceptions import ValidationError
from .schemas import Canonical, gwas_plot_spec
from .species import Species, resolve_species
from .utils import normalize_chrom, normalize_chrom_series
from .validation import check

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
        ValidationError: If neither species nor custom_order is provided, or
            the species is unknown or has no built-in chromosome order.
    """
    if custom_order is not None:
        return custom_order
    record = resolve_species(species)
    if record is None:
        raise ValidationError(
            "No chromosome order: pass a species, or custom_chrom_order in "
            "GenomeWideConfig"
        )
    if not record.chromosomes:
        raise ValidationError(
            f"No built-in chromosome order for species {record.key!r}; "
            f"pass custom_chrom_order in GenomeWideConfig."
        )
    return list(record.chromosomes)


def get_chromosome_colors(
    n_chromosomes: int, palette: Sequence[str] | None = None
) -> list[str]:
    """Get perceptually distinct colors for chromosomes.

    Uses a colorcet glasbey palette for good visual separation with
    saturated colors, unless the caller names one.

    Args:
        n_chromosomes: Number of chromosomes to color.
        palette: Colours to cycle through instead of the default.

    Returns:
        List of hex color strings.
    """
    if palette is None:
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
        gap: int = CHROMOSOME_GAP,
        palette: Sequence[str] | None = None,
    ) -> "GenomeLayout":
        """Lay out the genome axis once for every frame that shares it.

        Args:
            frames: Frames already filtered to plottable rows.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            order: Display order of the known chromosomes. Chromosomes the
                frames carry but this order omits are appended sorted.
            gap: Base pairs between one chromosome's end and the next's start.
            palette: Chromosome colours, or None for the default palette.

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
                cumulative += int(max_by_chrom[chrom]) + gap

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
            colors=dict(
                zip(full_order, get_chromosome_colors(len(full_order), palette))
            ),
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
    style: GenomeWideStyle = GenomeWideStyle(),
) -> list[PreparedManhattan]:
    """Validate each frame against ``config`` and lay them out on one genome.

    The boundary for the genome-wide families: every frame is checked for
    the chromosome, position and p-value columns the config names (and
    ``rs_col`` when given) before any of them is laid out, and projected onto
    the canonical columns with its chromosome names normalised.

    Args:
        dfs: GWAS results DataFrames, in panel order.
        config: Column names and chromosome order.
        species: Species whose chromosome order lays the axis out, unless
            ``config.custom_chrom_order`` overrides it.
        rs_col: SNP id column to require as well, or None.
        style: Supplies the chromosome gap and palette of the layout.

    Raises:
        ValidationError: If a frame is empty or lacks a named column.
    """
    normalized = []
    for df in dfs:
        check(
            df,
            gwas_plot_spec(
                config.pos_col, config.p_col, rs_col, chrom_col=config.chrom_col
            ),
        )
        roles = {
            Canonical.CHROM: normalize_chrom_series(df[config.chrom_col]),
            Canonical.POS: df[config.pos_col],
            Canonical.P: df[config.p_col],
        }
        if rs_col is not None:
            roles[Canonical.RS] = df[rs_col]
        normalized.append(pd.DataFrame(roles))
    custom_order = config.custom_chrom_order
    return prepare_manhattan_frames(
        normalized,
        species=species,
        custom_order=None
        if custom_order is None
        else [normalize_chrom(chrom) for chrom in custom_order],
        gap=style.chrom_gap,
        palette=style.palette,
    )


def prepare_manhattan_frames(
    dfs: Sequence[pd.DataFrame],
    *,
    species: str | Species | None = None,
    custom_order: list[str] | None = None,
    gap: int = CHROMOSOME_GAP,
    palette: Sequence[str] | None = None,
) -> list[PreparedManhattan]:
    """Lay out canonical GWAS frames against one shared genome layout.

    Every returned value carries the same :class:`GenomeLayout`, so a given
    ``(chrom, pos)`` lands at the same x in all of them. Pass a one-element
    list for a single panel. The frames carry the canonical ``chr``, ``pos``
    and ``p_value`` columns with normalised chromosome names, as
    :func:`prepare_genomewide_frames` projects them after validation.

    Args:
        dfs: Canonical GWAS frames, in panel order.
        species: Species for chromosome ordering.
        custom_order: Custom chromosome order.
        gap: Base pairs between one chromosome's end and the next's start.
        palette: Chromosome colours, or None for the default palette.

    Returns:
        One prepared value per input, in the same order, each carrying the
        columns ``_chrom_str``, ``_chrom_idx``, ``_cumulative_pos``,
        ``neglog10p`` and ``_color``.

    Raises:
        ValidationError: If no p-value in a frame survives, or if neither
            species nor custom_order names a chromosome order.
    """
    chrom_col, pos_col, p_col = Canonical.CHROM, Canonical.POS, Canonical.P
    order = get_chromosome_order(species, custom_order)
    filtered = [
        prepare_pvalue_data(
            df,
            p_col,
            "genome-wide",
            on_empty=ALL_PVALUES_INVALID.format(p_col=p_col),
        )
        for df in dfs
    ]
    layout = GenomeLayout.from_frames(
        filtered,
        chrom_col=chrom_col,
        pos_col=pos_col,
        order=order,
        gap=gap,
        palette=palette,
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
    palette: Sequence[str] | None = None,
) -> PreparedManhattan:
    """Prepare DataFrame for categorical Manhattan plot (PheWAS-style).

    Args:
        df: Results DataFrame with categories and p-values.
        category_col: Column name for category.
        p_col: Column name for p-value.
        category_order: Custom category order. Observed groups omitted from it
            are appended, sorted by label. Missing values form "Uncategorised".
        palette: Category colours, or None for the default palette.

    Returns:
        The frame with ``_cat_str``, ``_cat_idx``, ``_x_pos``, ``neglog10p``
        and ``_color`` columns, and its :class:`CategoryLayout`.
    """
    # Validate required columns
    if category_col not in df.columns:
        raise ValidationError(f"Column '{category_col}' not found in DataFrame")
    if p_col not in df.columns:
        raise ValidationError(f"Column '{p_col}' not found in DataFrame")

    result = prepare_pvalue_data(
        df, p_col, "genome-wide", on_empty=ALL_PVALUES_INVALID.format(p_col=p_col)
    )

    result["_cat_str"] = (
        result[category_col]
        .astype(object)
        .map(lambda value: "Uncategorised" if pd.isna(value) else str(value))
    )
    observed = set(result["_cat_str"])
    # An explicit order sets priority, but never hides observed groups.
    ordered = list(dict.fromkeys(str(value) for value in category_order or ()))
    category_order = ordered + sorted(observed - set(ordered))
    cat_to_idx = {cat: i for i, cat in enumerate(category_order)}
    result["_cat_idx"] = result["_cat_str"].map(cat_to_idx)

    # Use category index as x position (with jitter for multiple points per category)
    rng = np.random.default_rng(42)  # Local RNG for reproducible jitter
    result["_x_pos"] = result["_cat_idx"] + rng.uniform(-0.3, 0.3, size=len(result))

    layout = CategoryLayout(
        order=tuple(category_order),
        colors=dict(
            zip(category_order, get_chromosome_colors(len(category_order), palette))
        ),
    )
    result["_color"] = result["_cat_str"].map(layout.colors)
    return PreparedManhattan(result, layout)

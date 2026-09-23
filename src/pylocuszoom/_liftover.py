"""Coordinate-liftover port and the pure liftover transforms.

The liftover math is a pure per-position loop; its only external dependency is
a ``CoordinateLifter``, which is pyliftover's ``convert_coordinate`` contract.
``pyliftover.LiftOver`` satisfies it directly; tests use ``InMemoryLifter``, so
the coordinate math is exercised without a chain file or network access.

pyliftover positions are 0-based on both sides. Every position column in this
package is 1-based, so each query subtracts one and each hit adds it back.
"""

import os
import warnings
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from typing import List, Optional, Protocol, Tuple, Union, runtime_checkable

import pandas as pd

from .genome_build import GenomeBuild, resolve_build, ucsc_chrom
from .logging import logger
from .schemas import Canonical


@runtime_checkable
class CoordinateLifter(Protocol):
    """Maps a genomic coordinate to lifted coordinates in another build.

    This is pyliftover's ``LiftOver.convert_coordinate``, so a ``LiftOver``
    instance is a lifter as it stands.
    """

    def convert_coordinate(self, chrom: str, pos: int, /) -> Optional[List[Tuple]]:
        """Return ``[(chrom, pos, strand, score), ...]`` for a 0-based position.

        ``[]`` when the position does not lift; None when the chain does not
        know the chromosome at all.
        """
        ...


class InMemoryLifter:
    """Test adapter serving lifted coordinates from a dict.

    Keys are ``(chain chromosome, 0-based position)``. A value is a 0-based
    target position on the same chromosome, a ``(chrom, pos)`` tuple for a hit
    on another chromosome, or a list of positions for an ambiguous multi-mapping.
    A chromosome with no key returns None, as pyliftover does.
    """

    def __init__(self, mapping: dict) -> None:
        self._mapping = mapping
        self._chroms = {chrom for chrom, _ in mapping}

    def convert_coordinate(self, chrom: str, pos: int, /) -> Optional[List[Tuple]]:
        if chrom not in self._chroms:
            return None
        target = self._mapping.get((chrom, pos))
        if target is None:
            return []
        if isinstance(target, list):
            return [(chrom, p, "+", 0) for p in target]
        if isinstance(target, tuple):
            return [(target[0], target[1], "+", 0)]
        return [(chrom, target, "+", 0)]


@lru_cache(maxsize=4)
def load_chain(chain_path: Union[str, os.PathLike]) -> CoordinateLifter:
    """Load a UCSC chain file as a pyliftover lifter, once per path.

    Args:
        chain_path: Chain file, plain or gzipped.
    """
    from pyliftover import LiftOver

    logger.info(f"Loading liftover chain {chain_path}")
    return LiftOver(str(chain_path))


class _Outcome(Enum):
    LIFTED = "lifted"
    UNMAPPED = "unmapped"
    MULTI_MAPPED = "multi_mapped"
    CROSS_CHROM = "cross_chrom"


def _lift_one(
    lifter: CoordinateLifter, chrom: str, pos: int
) -> Tuple[_Outcome, Optional[int]]:
    """Lift one 1-based position, keeping only a single same-chromosome hit."""
    hits = lifter.convert_coordinate(chrom, int(pos) - 1)
    if not hits:
        return _Outcome.UNMAPPED, None
    if len(hits) > 1:
        return _Outcome.MULTI_MAPPED, None
    if hits[0][0] != chrom:
        return _Outcome.CROSS_CHROM, None
    return _Outcome.LIFTED, int(hits[0][1]) + 1


def liftover_positions(
    recomb_df: pd.DataFrame,
    lifter: CoordinateLifter,
    chrom: Optional[int] = None,
    build: Union[str, GenomeBuild, None] = None,
) -> pd.DataFrame:
    """Liftover every position in ``recomb_df`` through ``lifter``.

    Positions that do not lift to exactly one locus on the same chromosome are
    dropped; the result is sorted by position. Requires either a ``chr``
    column or the ``chrom`` argument.

    Raises:
        ValueError: If neither a ``chr`` column nor ``chrom`` is provided.
    """
    if "chr" in recomb_df.columns:
        chroms = recomb_df["chr"]
    elif chrom is not None:
        chroms = [chrom] * len(recomb_df)
    else:
        raise ValueError("Either 'chr' column or chrom parameter required")

    record = resolve_build(build)
    new_positions = [
        _lift_one(lifter, ucsc_chrom(chr_val, record), pos)[1]
        for chr_val, pos in zip(chroms, recomb_df["pos"])
    ]
    keep = [new_pos is not None for new_pos in new_positions]
    result_df = recomb_df[keep].copy()
    result_df["pos"] = [new_pos for new_pos in new_positions if new_pos is not None]

    unmapped = len(recomb_df) - len(result_df)
    if unmapped > 0:
        logger.debug(f"Dropped {unmapped} positions that failed to liftover")

    return result_df.sort_values("pos").reset_index(drop=True)


@dataclass(frozen=True, eq=False)
class RegionLiftResult:
    """Outcome of lifting one region's SNP positions to another genome build.

    Attributes:
        lifted_df: The input rows that lifted, in input order, with the
            position column replaced by the target-build position. Every other
            column, p-values and LD included, travels with its row unchanged.
        lead_pos: Lifted lead position, or None when no lead was asked for or
            it did not lift cleanly.
        source_lead_pos: The lead position the caller asked to lift.
        start: Smallest lifted position, or None when nothing lifted.
        end: Largest lifted position, or None when nothing lifted.
        n_unmapped: Rows the chain has no target for.
        n_multimapped: Rows with more than one target locus.
        n_cross_chrom: Rows whose only target is on another chromosome.
        is_collinear: False when the lifted SNPs change order between builds,
            meaning the region is locally rearranged and a left-to-right plot
            misrepresents it.
    """

    lifted_df: pd.DataFrame
    lead_pos: Optional[int]
    source_lead_pos: Optional[int]
    start: Optional[int]
    end: Optional[int]
    n_unmapped: int
    n_multimapped: int
    n_cross_chrom: int
    is_collinear: bool

    @property
    def n_lifted(self) -> int:
        """Rows that lifted to exactly one locus on the same chromosome."""
        return len(self.lifted_df)

    @property
    def n_dropped(self) -> int:
        """Rows dropped for any reason."""
        return self.n_unmapped + self.n_multimapped + self.n_cross_chrom

    @property
    def n_input(self) -> int:
        """Rows the caller asked to lift."""
        return self.n_lifted + self.n_dropped


def liftover_region(
    region_df: pd.DataFrame,
    *,
    chrom: Union[int, str],
    lifter: CoordinateLifter,
    pos_col: str = Canonical.POS,
    lead_pos: Optional[int] = None,
    build: Union[str, GenomeBuild, None] = None,
) -> RegionLiftResult:
    """Lift one region's SNP positions to another genome build.

    A SNP is kept only when it lifts to exactly one locus on the same
    chromosome; unmapped, multi-mapped and cross-chromosome SNPs are dropped
    and counted by reason. Only positions change.

    Args:
        region_df: GWAS rows for one chromosome, in the source build.
        chrom: Chromosome of the region, bare or ``chr``-prefixed. The chain
            is queried with the UCSC name, so ``build``'s renames apply.
        lifter: Coordinate lifter exposing ``convert_coordinate``, such as
            ``pyliftover.LiftOver``.
        pos_col: 1-based position column to lift.
        lead_pos: 1-based lead-SNP position to lift alongside the region.
        build: Build whose UCSC chromosome renames apply, such as canine
            PLINK code 39 queried as ``chrX``; a name or a record. None, or a
            build the package does not know, applies none.

    Returns:
        RegionLiftResult with the lifted rows and the drop counts.

    Example:
        >>> from pyliftover import LiftOver
        >>> lift = liftover_region(
        ...     region_df, chrom=12, lifter=LiftOver("canFam3ToCanFam4.over.chain.gz")
        ... )
        >>> lift.n_lifted, lift.n_dropped
    """
    name = ucsc_chrom(chrom, resolve_build(build))
    if len(region_df) and lifter.convert_coordinate(name, 0) is None:
        warnings.warn(
            f"{name} is unknown to the liftover chain, so no SNP in this region "
            "can lift; check chromosome naming against the chain",
            stacklevel=2,
        )

    outcomes = [_lift_one(lifter, name, pos) for pos in region_df[pos_col]]
    keep = [outcome is _Outcome.LIFTED for outcome, _ in outcomes]
    lifted = region_df[keep].copy()
    lifted[pos_col] = pd.Series(
        [new_pos for _, new_pos in outcomes if new_pos is not None],
        index=lifted.index,
        dtype="int64",
    )

    counts = Counter(outcome for outcome, _ in outcomes)

    old_kept = region_df.loc[keep, pos_col].to_numpy()
    new_kept = lifted[pos_col].to_numpy()
    in_old_order = new_kept[old_kept.argsort(kind="stable")]
    is_collinear = bool((in_old_order[1:] >= in_old_order[:-1]).all())

    return RegionLiftResult(
        lifted_df=lifted,
        lead_pos=(
            _lift_one(lifter, name, lead_pos)[1] if lead_pos is not None else None
        ),
        source_lead_pos=lead_pos,
        start=int(new_kept.min()) if len(new_kept) else None,
        end=int(new_kept.max()) if len(new_kept) else None,
        n_unmapped=counts[_Outcome.UNMAPPED],
        n_multimapped=counts[_Outcome.MULTI_MAPPED],
        n_cross_chrom=counts[_Outcome.CROSS_CHROM],
        is_collinear=is_collinear,
    )

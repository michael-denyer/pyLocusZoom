"""Coordinate-liftover port and the pure liftover transforms.

The liftover math is a pure per-position loop; its only external dependency is
a ``CoordinateLifter``, which is pyliftover's ``convert_coordinate`` contract.
``pyliftover.LiftOver`` satisfies it directly; tests use ``InMemoryLifter``, so
the coordinate math is exercised without a chain file or network access.

pyliftover positions are 0-based on both sides. Every position column in this
package is 1-based, so each query subtracts one and each hit adds it back.
"""

import gzip
import os
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Protocol, Tuple, Union, runtime_checkable

import pandas as pd

from ._http import download_file
from .exceptions import DataDownloadError, OptionalDependencyMissing, ValidationError
from .genome_build import GenomeBuild, resolve_build, ucsc_chrom
from .logging import logger
from .schemas import Canonical
from .utils import _platform_cache_base, filter_by_region


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

    The only place a chain file becomes a lifter, whether the caller named it
    or the package downloaded it.

    Args:
        chain_path: Chain file, plain or gzipped.

    Raises:
        OptionalDependencyMissing: If pyliftover is not installed.
        ValidationError: If the chain is missing, unreadable or has no mappings.
    """
    try:
        from pyliftover import LiftOver
    except ImportError as e:
        raise OptionalDependencyMissing(
            "pyliftover is required for liftover. Install it with: pip install pyliftover"
        ) from e

    logger.info(f"Loading liftover chain {chain_path}")
    path = Path(chain_path)
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        # pyliftover raises bare Exception for malformed headers and blocks.
        # Translate at this adapter, and close the stream even when parsing fails.
        with opener(path, "rb") as stream:
            lifter = LiftOver(stream)
    except Exception as e:
        raise ValidationError(f"Liftover chain {path} is unreadable: {e}") from e
    if not lifter.chain_file.chains:
        raise ValidationError(f"Liftover chain {path} contains no mappings")
    return lifter


def chain_lifter(source: GenomeBuild, target: GenomeBuild) -> CoordinateLifter:
    """Return the lifter for a registered chain, downloading it on first use.

    Chains are cached under the platform cache's ``liftover`` leaf, beside the
    recombination maps rather than inside them. A cached chain that no longer
    parses is downloaded once more before giving up.

    Args:
        source: Build the coordinates are in.
        target: Build to lift them to.

    Raises:
        ValidationError: If no chain from ``source`` to ``target`` is registered.
        DataDownloadError: If the chain cannot be downloaded, or is unreadable
            after a fresh download.
        OptionalDependencyMissing: If pyliftover is not installed.
    """
    url = source.chain_url(target)
    if url is None:
        raise ValidationError(
            f"No liftover chain is registered from {source.assembly_name} "
            f"to {target.assembly_name}"
        )
    path = _platform_cache_base() / "liftover" / url.rsplit("/", 1)[-1]
    if not path.exists():
        _download_chain(url, path)
    try:
        return load_chain(path)
    except ValidationError as e:
        logger.warning(f"Liftover chain {path} is unreadable ({e}); refetching")
    _download_chain(url, path)
    try:
        return load_chain(path)
    except ValidationError as e:
        try:
            path.unlink(missing_ok=True)
        except OSError as cleanup_error:
            raise DataDownloadError(
                f"Liftover chain {path} is unreadable and could not be removed: {cleanup_error}"
            ) from e
        raise DataDownloadError(f"Liftover chain {path} is unreadable: {e}") from e


def _download_chain(url: str, path: Path) -> None:
    """Download one chain file into the cache."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Downloading liftover chain {url}")
        download_file(url, path, desc="Liftover chain")
    except OSError as e:
        raise DataDownloadError(f"Could not write liftover chain {path}: {e}") from e


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
        chain_has_chrom: False when the chain does not know the region's
            chromosome at all, so nothing could lift; usually a naming
            mismatch between the frame and the chain.
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
    chain_has_chrom: bool = True

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
        chain_has_chrom=lifter.convert_coordinate(name, 0) is not None,
    )


def describe_drops(result: RegionLiftResult, chrom: Union[int, str]) -> str:
    """Say why rows did not lift, naming a chromosome the chain lacks."""
    if not result.chain_has_chrom:
        return (
            f"chr{chrom} is unknown to the liftover chain; check chromosome "
            "naming against the chain"
        )
    return (
        f"{result.n_unmapped} unmapped, {result.n_multimapped} multi-mapped, "
        f"{result.n_cross_chrom} on another chromosome"
    )


@dataclass(frozen=True)
class LiftedWindow:
    """Regional frames lifted to the target build, and the window around them.

    Attributes:
        frames: Each input frame's in-region rows that lifted, in input order.
        start: Window start in the target build, keeping the requested margin
            before the first lifted SNP of any frame.
        end: Window end in the target build, keeping the requested margin
            after the last lifted SNP of any frame.
        lead_positions: Each requested lead lifted, or None where there was
            none or it did not lift.
        notes: Sentences the caller should pass on to the user: a lead that
            did not lift, or a region rearranged between builds.
    """

    frames: List[pd.DataFrame]
    start: int
    end: int
    lead_positions: List[Optional[int]]
    notes: Tuple[str, ...]


def lift_window(
    frames: List[pd.DataFrame],
    *,
    chrom: Union[int, str],
    start: int,
    end: int,
    chrom_col: Optional[str],
    pos_col: str,
    lead_positions: List[Optional[int]],
    lifter: CoordinateLifter,
    build: Union[str, GenomeBuild, None],
) -> LiftedWindow:
    """Lift one region of each frame to another build, with the window around them.

    Every regional entry point lifts through this, so ``plot()`` and
    ``plot_stacked()`` cannot drift. The requested bounds need not lift
    themselves, so the window keeps the requested margins around the
    outermost lifted SNPs; with several frames it spans all of them.

    Args:
        frames: Source-build frames, one per panel.
        chrom: Chromosome of the region.
        start: Source-build region start.
        end: Source-build region end.
        chrom_col: Chromosome column in every frame; None selects the region
            by position only, for frames already scoped to one chromosome.
        pos_col: 1-based position column in every frame.
        lead_positions: One source-build lead per frame, or None for none.
        lifter: Lifter from the source build to the target build.
        build: Target build, used for its chromosome renames and in messages.

    Raises:
        ValidationError: If no SNP of a frame's region lifts.
    """
    target = resolve_build(build)
    target_name = target.assembly_name if target else (build or "the target build")
    where = f"chr{chrom}:{start}-{end}"
    lifted_frames, leads, starts, ends, notes = [], [], [], [], []
    for frame, lead_pos in zip(frames, lead_positions):
        selected = filter_by_region(
            frame, region=(chrom, start, end), chrom_col=chrom_col, pos_col=pos_col
        )
        lift = liftover_region(
            selected,
            chrom=chrom,
            lifter=lifter,
            pos_col=pos_col,
            lead_pos=lead_pos,
            build=target,
        )
        if lift.lifted_df.empty:
            raise ValidationError(
                f"No SNP in {where} lifted to {target_name}: "
                f"{describe_drops(lift, chrom)}"
            )
        if lift.n_dropped:
            logger.info(
                "Liftover dropped {}/{} SNPs in {} ({})",
                lift.n_dropped,
                lift.n_input,
                where,
                describe_drops(lift, chrom),
            )
        if not lift.is_collinear:
            notes.append(
                f"{where} is rearranged between builds; the regional plot's "
                "left-to-right order may misrepresent it"
            )
        if lead_pos is not None and lift.lead_pos is None:
            notes.append(
                f"Lead SNP at chr{chrom}:{lead_pos} did not lift to "
                f"{target_name}; the lead is auto-detected instead"
            )
        source_pos = selected.loc[lift.lifted_df.index, pos_col]
        window_start = max(1, lift.start - int(source_pos.min() - start))
        starts.append(window_start)
        ends.append(max(lift.end + int(end - source_pos.max()), window_start + 1))
        lifted_frames.append(lift.lifted_df)
        leads.append(lift.lead_pos)
    return LiftedWindow(
        frames=lifted_frames,
        start=min(starts),
        end=max(ends),
        lead_positions=leads,
        notes=tuple(dict.fromkeys(notes)),
    )

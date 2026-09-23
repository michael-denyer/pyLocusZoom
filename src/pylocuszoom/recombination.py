"""Recombination rate overlay and data management.

Provides:
- Recombination rate overlay for regional plots
- Download and loading of species-specific recombination maps
- Liftover of the maps through the chains registered on their GenomeBuild
"""

import io
import os
import re
import shutil
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Optional

import pandas as pd

from ._http import download_file
from ._liftover import CoordinateLifter, chain_lifter, describe_drops, liftover_region
from .exceptions import DataDownloadError, RecombinationMapNotFound, ValidationError
from .genome_build import GENOME_BUILDS, assembly_token, resolve_build
from .logging import logger
from .species import Species, resolve_species
from .utils import _platform_cache_base, filter_by_region

CANINE_MAP_FILENAMES = frozenset(f"chr{chrom}_recomb.tsv" for chrom in range(1, 39))

# Data sources by species
CANINE_RECOMB_URL = (
    "https://github.com/cflerin/dog_recombination/raw/master/dog_genetic_maps.tar.gz"
)


@dataclass(frozen=True)
class RecombSource:
    """Where one species' built-in recombination maps come from.

    Everything that varies between species lives here, so
    ``download_recombination_maps`` has no species-specific branch and a second
    species is one more row rather than a second downloader.

    Attributes:
        species: Canonical Species.key this source serves.
        url: Archive holding the published map set.
        archive_glob: Glob selecting regular map members inside the
            archive. Deliberately loose; chrom_pattern is what decides
            whether a matched file really is a map.
        chrom_pattern: Regex anchored at the start of a matched file's stem
            whose first group is the chromosome name. A matched file this
            does not parse is a DataDownloadError, not a file to skip: the
            old four-stage string peel turned an unrecognised name into a
            silent chr_recomb.tsv.
        filenames: The complete map set. A directory holding exactly these is
            a cache hit; anything else is downloaded again.
        native_build: Key of the GenomeBuild the published maps are in. Its
            ``liftover_chains`` name the builds the maps can be lifted to;
            other non-native builds are unavailable.
    """

    species: str
    url: str
    archive_glob: str
    chrom_pattern: str
    filenames: frozenset[str]
    native_build: str


CANINE_SOURCE = RecombSource(
    species="canine",
    url=CANINE_RECOMB_URL,
    archive_glob="chr*.txt",
    chrom_pattern=r"chr(\d+|X|Y|MT)(?:_|$)",
    filenames=CANINE_MAP_FILENAMES,
    native_build="canfam3",
)

# A species absent from here has no built-in maps and supplies its own.
RECOMB_SOURCES: dict[str, RecombSource] = {CANINE_SOURCE.species: CANINE_SOURCE}


_CANONICAL_RECOMB_HEADER = "chr\tpos\trate\tcM\n"
_KNOWN_HEADER_TOKENS = frozenset(
    {"chr", "chrom", "chromosome", "pos", "position", "bp"}
)


def ensure_recomb_header(content: str, source_name: str) -> str:
    r"""Return content with a canonical header row, prepending one if absent.

    A numeric first token means the first row is data, so a
    ``chr\tpos\trate\tcM`` header is prepended. A non-numeric first token must be
    one of the known header names; anything else (e.g. ``<html>`` from a
    corrupted mirror or an HTTP error body) is rejected rather than silently
    treated as a header.

    Args:
        content: Raw text of a recombination map file.
        source_name: Source file name, used in the error message.

    Returns:
        The content, with a header row prepended if it was missing.

    Raises:
        DataDownloadError: If the first token is non-numeric and not a known
            header.
    """
    lines = content.strip().split("\n")
    first_token = lines[0].split()[0] if lines[0].split() else ""
    normalised_token = first_token.lstrip("#").lower()
    try:
        float(first_token)
        has_header = False
    except ValueError:
        if normalised_token not in _KNOWN_HEADER_TOKENS:
            raise DataDownloadError(
                f"Unrecognised first token {first_token!r} in recombination "
                f"map {source_name}; refusing to treat as header. The "
                f"downloaded archive may be corrupted."
            )
        has_header = True
    if not has_header:
        content = _CANONICAL_RECOMB_HEADER + content
    return content


def get_default_data_dir() -> Path:
    """Get default directory for recombination map data.

    Returns the shared platform cache root (see ``_platform_cache_base``) with
    the ``recombination_maps`` leaf appended:
    - macOS/Linux: ~/.cache/pylocuszoom/recombination_maps (or $XDG_CACHE_HOME)
    - Windows: %LOCALAPPDATA%/pylocuszoom/recombination_maps
    - Databricks: /dbfs/FileStore/reference_data/recombination_maps
    """
    return _platform_cache_base() / "recombination_maps"


def _resolve_map_dir(data_dir: str | Path | None) -> Path:
    """Return the directory recombination maps live in, defaulting to the cache."""
    return get_default_data_dir() if data_dir is None else Path(data_dir)


def _has_complete_maps(path: Path, source: RecombSource) -> bool:
    """Return whether path contains exactly this source's map set."""
    if not path.exists():
        return False
    present = {map_path.name for map_path in path.glob("chr*_recomb.tsv")}
    return present == source.filenames


def _holds_only_maps(path: Path, source: RecombSource) -> bool:
    """Return whether replacing path wholesale would discard only this source's maps."""
    if not path.exists():
        return True
    return path.is_dir() and all(
        entry.is_file() and entry.name in source.filenames for entry in path.iterdir()
    )


def _publish_map_generation(
    staging_dir: Path, output_path: Path, source: RecombSource
) -> Path:
    """Install a complete map set without the directory ever going missing.

    An absent target receives the staging directory in one rename. A target
    that exists (a set being refreshed, a damaged set, or one a concurrent
    writer installed first) keeps its directory while each map file is
    swapped in with one ``os.replace``, so a reader finds the old file or the
    new one, never a gap. Nothing is moved aside, so nothing can be left
    behind, and a writer that loses the race to another still succeeds. The
    map set for a source URL never changes, so interleaved writers converge
    on the same files.
    """
    if not _has_complete_maps(staging_dir, source):
        raise DataDownloadError(
            f"Downloaded recombination archive does not contain the complete "
            f"{source.species} map set ({len(source.filenames)} chromosome files)"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.is_symlink():
        # Older releases published behind a symlink; its target is not ours.
        output_path.unlink()
    try:
        os.rename(staging_dir, output_path)
        return output_path
    except OSError:
        if not output_path.is_dir():
            raise
    for name in sorted(source.filenames):
        os.replace(staging_dir / name, output_path / name)
    for stray in output_path.glob("chr*_recomb.tsv"):
        if stray.name not in source.filenames:
            stray.unlink(missing_ok=True)
    return output_path


def _stage_archive(tar_path: Path, source: RecombSource, staging: Path) -> None:
    """Read regular map members into canonical files without extracting paths."""
    pattern = re.compile(source.chrom_pattern, re.IGNORECASE)
    seen = set()
    try:
        with tarfile.open(tar_path, "r:gz") as archive:
            for member in archive:
                name = PurePosixPath(member.name)
                if name.is_absolute() or ".." in name.parts:
                    raise DataDownloadError(
                        f"Archive contains unsafe member {member.name!r}"
                    )
                if member.isdir():
                    continue
                if not member.isfile():
                    raise DataDownloadError(
                        f"Archive member {member.name!r} is not a regular file"
                    )
                if not name.match(source.archive_glob):
                    continue
                match = pattern.match(name.stem)
                if match is None:
                    raise DataDownloadError(
                        f"Recombination map {name.name!r} does not name a chromosome; "
                        f"expected a stem matching {source.chrom_pattern!r}"
                    )
                filename = f"chr{match.group(1)}_recomb.tsv"
                if filename not in source.filenames:
                    raise DataDownloadError(f"Unexpected chromosome map {name.name!r}")
                if filename in seen:
                    raise DataDownloadError(f"Duplicate chromosome map {filename!r}")
                seen.add(filename)
                with (
                    io.TextIOWrapper(
                        archive.extractfile(member), encoding="utf-8"
                    ) as body,
                    (staging / filename).open("w", encoding="utf-8") as output,
                ):
                    output.write(ensure_recomb_header(body.readline(), member.name))
                    shutil.copyfileobj(body, output)
    except (tarfile.TarError, UnicodeDecodeError) as e:
        raise DataDownloadError(
            f"Downloaded recombination archive from {source.url} is not a valid tar.gz map set: {e}"
        ) from e
    if not seen:
        raise DataDownloadError(
            f"Could not find chromosome map files matching "
            f"{source.archive_glob!r} in archive from {source.url}"
        )


def download_recombination_maps(source: RecombSource, output_path: Path) -> Path:
    """Download, extract and publish one source's complete map set.

    Unconditional: the caller owns the cache-hit decision. Everything is
    written into a temporary directory and promoted with one rename, so a
    failure part-way through cannot leave a partial set behind that a later
    cache check would accept.

    Args:
        source: Source to download.
        output_path: Directory the complete map set is published to.

    Returns:
        ``output_path``.

    Raises:
        DataDownloadError: If the download fails, the archive is corrupt,
            incomplete, or not a recombination map set, or the maps cannot be
            written.
    """
    logger.info(f"Downloading {source.species} recombination maps...")
    logger.debug(f"Source: {source.url}")

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output_path.parent) as tmpdir:
            tmp = Path(tmpdir)
            archive = tmp / "maps.tar.gz"
            download_file(source.url, archive, desc="Recombination maps")
            logger.debug(f"Downloaded {archive.stat().st_size / 1024:.1f} KB")

            staging = tmp / "_staging"
            staging.mkdir()
            _stage_archive(archive, source, staging)
            _publish_map_generation(staging, output_path, source)
    except OSError as e:
        raise DataDownloadError(
            f"Could not write recombination maps to {output_path}: {e}"
        ) from e

    logger.info(f"Recombination maps saved to: {output_path}")
    return output_path


def download_canine_recombination_maps(
    output_dir: Optional[str] = None,
    force: bool = False,
) -> Path:
    """Download canine recombination rate maps from Campbell et al. 2016.

    Downloads from: https://github.com/cflerin/dog_recombination

    Data is in CanFam3.1 coordinates with columns:
    - chr: Chromosome number
    - pos: Physical position (bp)
    - rate: Recombination rate (cM/Mb)
    - cM: Cumulative genetic distance (centiMorgans)

    Args:
        output_dir: Directory to save maps. Uses platform cache if None. It
            must be new, empty or hold only a previous canine map set, because
            publishing replaces the whole directory.
        force: Re-download even if files exist.

    Returns:
        Path to the directory containing recombination map files.

    Raises:
        DataDownloadError: If the download fails or the archive is corrupt,
            incomplete, or not a recombination map set.
        ValidationError: If output_dir holds anything besides canine maps.
    """
    output_path = _resolve_map_dir(output_dir)
    if not force and _has_complete_maps(output_path, CANINE_SOURCE):
        return output_path
    if output_dir is not None and not _holds_only_maps(output_path, CANINE_SOURCE):
        raise ValidationError(
            f"output_dir {output_path} holds files other than the canine "
            "recombination maps, and publishing the map set would replace it. "
            "Pass a new or empty directory."
        )
    return download_recombination_maps(CANINE_SOURCE, output_path)


def load_recombination_map(
    chrom: int,
    species: str | Species | None = "canine",
    data_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Load recombination map for a specific chromosome.

    Args:
        chrom: Chromosome number (1-38 for canine, 1-18 for feline) or 'X'.
        species: Species name, alias or record, or None for a caller
            supplying its own maps.
        data_dir: Directory containing recombination maps.

    Returns:
        DataFrame with columns: pos, rate, cM.

    Raises:
        RecombinationMapNotFound: If the map file is not there.
        ValidationError: If the species is not one this package knows.
    """
    record = resolve_species(species)
    data_path = _resolve_map_dir(data_dir)
    chrom_str = str(chrom).replace("chr", "")
    map_file = data_path / f"chr{chrom_str}_recomb.tsv"

    if not map_file.exists():
        key = record.key if record else None
        remedy = (
            f"Run ensure_recomb_maps(species={key!r}) first to download them."
            if key in RECOMB_SOURCES
            else f"There are no built-in recombination maps for {key!r}; "
            f"pass data_dir pointing at your own chr{{N}}_recomb.tsv files."
        )
        raise RecombinationMapNotFound(
            f"Recombination map not found: {map_file}\n{remedy}"
        )

    df = pd.read_csv(map_file, sep="\t")

    for col in ("pos", "rate", "cM"):
        if col not in df.columns:
            continue
        original = df[col]
        df[col] = pd.to_numeric(df[col], errors="coerce")
        dropped = original[df[col].isna() & original.notna()]
        if not dropped.empty:
            logger.warning(
                f"Recombination map chr{chrom_str}: {len(dropped)} non-numeric "
                f"values in '{col}' column dropped. "
                f"Sample values: {dropped.head(3).tolist()}"
            )

    return df.dropna(subset=["pos", "rate"])


def get_recombination_rate_for_region(
    chrom: int,
    start: int,
    end: int,
    species: str | Species | None = "canine",
    data_dir: Optional[str] = None,
    genome_build: Optional[str] = None,
    lifter: Optional[CoordinateLifter] = None,
) -> pd.DataFrame:
    """Get recombination rate data for a genomic region, or raise why there is none.

    The one lookup behind the recombination overlay. Every reason there is
    no frame is a ``PyLocusZoomError`` subclass, so a caller that would rather
    draw without the overlay catches that one base class; the plotter does.

    Args:
        chrom: Chromosome number.
        start: Start position (bp).
        end: End position (bp).
        species: Species name, alias or record, or None for a caller
            supplying its own maps.
        data_dir: Read-only caller maps, already in the target build. None
            uses the managed maps, downloading them on first use.
        genome_build: Target genome build (e.g., "canfam4"). Managed maps
            use the chain their build registers when conversion is needed.
            Caller maps are already in this build.
        lifter: Lift the loaded maps, managed or caller, through this lifter
            instead. The maps must be in the lifter's source build; the
            registered chain is then not consulted, and ``genome_build``
            only supplies the UCSC chromosome names.

    Returns:
        DataFrame with pos and rate columns for the region.

    Raises:
        RecombinationMapNotFound: If the species has no built-in maps and no
            ``data_dir`` was given, or there is no map for the chromosome.
        DataDownloadError: If the managed maps or the liftover chain cannot
            be downloaded, read or written.
        OptionalDependencyMissing: If liftover needs pyliftover and it is
            not installed.
        ValidationError: If the species is unknown, no chain reaches
            ``genome_build``, or the lift maps none of the chromosome's
            positions.

    Note:
        Built-in canine recombination maps are in CanFam3.1 coordinates.
        With no data_dir and genome_build="canfam4", positions are lifted over.
        This requires pyliftover: pip install pyliftover
    """
    record = resolve_species(species)
    map_dir = ensure_recomb_maps(species=record, data_dir=data_dir)
    if map_dir is None:
        raise RecombinationMapNotFound(
            f"There are no built-in recombination maps for "
            f"{record.key if record else None!r}; pass data_dir pointing at "
            "your own chr{N}_recomb.tsv files."
        )
    source = RECOMB_SOURCES[record.key] if data_dir is None and lifter is None else None
    target = resolve_build(genome_build)
    lift_build = target
    if (
        source is not None
        and genome_build
        and assembly_token(genome_build) != source.native_build
    ):
        native = GENOME_BUILDS[source.native_build]
        if target is None or native.chain_url(target) is None:
            raise ValidationError(
                f"Built-in {source.species} maps use {source.native_build}; "
                f"no liftover chain is available for {genome_build!r}. "
                "Supply data_dir with maps in the requested build."
            )
        lifter, lift_build = chain_lifter(native, target), native
    df = load_recombination_map(chrom, species=record, data_dir=map_dir)
    if lifter is not None:
        logger.debug(f"Lifting over recombination map for chr{chrom}")
        lift = liftover_region(
            df, chrom=chrom, lifter=lifter, pos_col="pos", build=lift_build
        )
        if lift.lifted_df.empty and not df.empty:
            raise ValidationError(
                f"Liftover mapped none of the {len(df)} positions in the "
                f"chr{chrom} recombination map: {describe_drops(lift, chrom)}"
            )
        df = lift.lifted_df.sort_values("pos").reset_index(drop=True)

    # Filter to region
    region_df = filter_by_region(
        df,
        region=(chrom, start, end),
        chrom_col=None,  # A recombination map has no chromosome column
        pos_col="pos",
    )

    return region_df[["pos", "rate"]]


def ensure_recomb_maps(
    species: str | Species | None = "canine",
    data_dir: Optional[str] = None,
) -> Optional[Path]:
    """Ensure recombination maps are available, downloading if needed.

    Args:
        species: Species name, alias or record, or None for a caller
            supplying its own maps.
        data_dir: Read-only caller directory. Only None selects the managed cache
            and permits downloads.

    Returns:
        Path to the recombination maps directory, or None if the species has
        no built-in map set and no caller directory was supplied.

    Raises:
        ValidationError: If the species is not one this package knows.
        DataDownloadError: If the species has maps and they could not be
            fetched or written.
    """
    record = resolve_species(species)
    if data_dir is not None:
        return Path(data_dir)
    source = RECOMB_SOURCES.get(record.key) if record else None
    if source is None:
        logger.debug(f"No built-in recombination maps for species: {species}")
        return None

    output_path = get_default_data_dir()

    if _has_complete_maps(output_path, source):
        logger.debug(f"Recombination maps already exist at {output_path}")
        return output_path

    return download_recombination_maps(source, output_path)

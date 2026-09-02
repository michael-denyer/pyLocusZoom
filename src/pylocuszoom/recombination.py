"""Recombination rate overlay and data management.

Provides:
- Recombination rate overlay for regional plots
- Download and loading of species-specific recombination maps
- Liftover support for CanFam3.1 to CanFam4 coordinate conversion
"""

import os
import re
import shutil
import tarfile
import tempfile
import uuid
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

import pandas as pd

from ._http import download_file
from ._liftover import PyLiftOverLifter, liftover_positions
from .exceptions import DataDownloadError, OptionalDependencyMissing
from .logging import logger
from .species import Species, resolve_species
from .utils import _platform_cache_base, assembly_token, filter_by_region

# Recombination overlay color
RECOMB_COLOR = "#7FCDFF"  # Light blue
CANINE_MAP_FILENAMES = frozenset(f"chr{chrom}_recomb.tsv" for chrom in range(1, 39))

# Data sources by species
CANINE_RECOMB_URL = (
    "https://github.com/cflerin/dog_recombination/raw/master/dog_genetic_maps.tar.gz"
)

# Liftover chain files
CANFAM3_TO_CANFAM4_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/gbdb/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz"


@dataclass(frozen=True)
class RecombSource:
    """Where one species' built-in recombination maps come from.

    Everything that varies between species lives here, so
    ``download_recombination_maps`` has no species-specific branch and a second
    species is one more row rather than a second downloader.

    Attributes:
        species: Canonical Species.key this source serves.
        url: Archive holding the published map set.
        archive_glob: Glob selecting the map files inside the extracted
            archive. Deliberately loose; chrom_pattern is what decides
            whether a matched file really is a map.
        chrom_pattern: Regex anchored at the start of a matched file's stem
            whose first group is the chromosome name. A matched file this
            does not parse is a DataDownloadError, not a file to skip: the
            old four-stage string peel turned an unrecognised name into a
            silent chr_recomb.tsv.
        filenames: The complete map set. A directory holding exactly these is
            a cache hit; anything else is downloaded again.
        native_build: Build the published maps are in.
        liftover_chains: Target build token -> chain URL, for the builds the
            maps can be lifted to. A build absent from here is served in
            ``native_build`` coordinates.
    """

    species: str
    url: str
    archive_glob: str
    chrom_pattern: str
    filenames: frozenset[str]
    native_build: str
    liftover_chains: dict[str, str]


CANINE_SOURCE = RecombSource(
    species="canine",
    url=CANINE_RECOMB_URL,
    archive_glob="chr*.txt",
    chrom_pattern=r"chr(\d+|X|Y|MT)(?:_|$)",
    filenames=CANINE_MAP_FILENAMES,
    native_build="canfam3",
    liftover_chains={"canfam4": CANFAM3_TO_CANFAM4_CHAIN_URL},
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


def get_chain_dir() -> Path:
    """Get the directory holding downloaded liftover chain files.

    Chains live beside the recombination maps rather than inside them, so a
    map set can be replaced wholesale without taking the chain with it.
    """
    return _platform_cache_base() / "liftover"


def get_chain_file_path() -> Path:
    """Get path to the CanFam3 to CanFam4 liftover chain file."""
    return get_chain_dir() / CANFAM3_TO_CANFAM4_CHAIN_URL.rsplit("/", 1)[-1]


def download_liftover_chain(force: bool = False) -> Path:
    """Download the CanFam3 to CanFam4 liftover chain file.

    Args:
        force: Re-download even if file exists.

    Returns:
        Path to the downloaded chain file.

    Raises:
        DataDownloadError: If the download fails.
    """
    chain_path = get_chain_file_path()

    if chain_path.exists() and not force:
        return chain_path

    chain_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading CanFam3 to CanFam4 liftover chain...")
    logger.debug(f"Source: {CANFAM3_TO_CANFAM4_CHAIN_URL}")

    download_file(
        CANFAM3_TO_CANFAM4_CHAIN_URL,
        chain_path,
        desc="Liftover chain",
    )

    logger.info(f"Chain file saved to: {chain_path}")
    return chain_path


def liftover_recombination_map(
    recomb_df: pd.DataFrame,
    from_build: str = "canfam3",
    to_build: str = "canfam4",
    chrom: Optional[int] = None,
) -> pd.DataFrame:
    """Liftover recombination map coordinates between genome builds.

    Args:
        recomb_df: DataFrame with 'pos' column (and optionally 'chr').
        from_build: Source genome build (default: canfam3).
        to_build: Target genome build (default: canfam4).
        chrom: Chromosome number (required if 'chr' not in recomb_df).

    Returns:
        DataFrame with lifted coordinates. Positions that fail to map are dropped.
    """
    try:
        import pyliftover  # noqa: F401
    except ImportError as e:
        raise OptionalDependencyMissing(
            "pyliftover is required for CanFam4 liftover. "
            "Install it with: pip install pyliftover"
        ) from e

    chain_path = download_liftover_chain()
    logger.debug(f"Lifting over coordinates from {from_build} to {to_build}")
    return liftover_positions(recomb_df, PyLiftOverLifter(chain_path), chrom)


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


def _discard_map_dir(path: Path) -> None:
    """Delete a replaced map directory, following the link older releases left."""
    if path.is_symlink():
        target = path.resolve()
        path.unlink()
        if target.is_dir():
            shutil.rmtree(target, ignore_errors=True)
    else:
        shutil.rmtree(path, ignore_errors=True)


def _publish_map_generation(
    staging_dir: Path, output_path: Path, source: RecombSource
) -> Path:
    """Swap a complete map set into place, keeping the old one until it lands."""
    if not _has_complete_maps(staging_dir, source):
        raise DataDownloadError(
            f"Downloaded recombination archive does not contain the complete "
            f"{source.species} map set ({len(source.filenames)} chromosome files)"
        )

    parent = output_path.parent
    parent.mkdir(parents=True, exist_ok=True)
    previous = parent / f".{output_path.name}.previous-{uuid.uuid4().hex}"
    replacing = output_path.is_symlink() or output_path.exists()
    if replacing:
        os.replace(output_path, previous)

    try:
        os.replace(staging_dir, output_path)
    except BaseException:
        if replacing:
            os.replace(previous, output_path)
        raise

    if replacing:
        _discard_map_dir(previous)
    return output_path


def _extract_archive(tar_path: Path, into: Path, url: str) -> Path:
    """Extract a downloaded tar.gz, refusing members that escape ``into``.

    Args:
        tar_path: The downloaded archive.
        into: Directory to extract beneath.
        url: Source URL, named in the error if the archive is not a tar.gz.

    Returns:
        The directory the archive was extracted into.

    Raises:
        DataDownloadError: If the file is not a readable tar.gz.
    """
    root = into.resolve()
    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            safe_members = []
            for member in tar.getmembers():
                try:
                    (into / member.name).resolve().relative_to(root)
                    safe_members.append(member)
                except (ValueError, OSError):
                    logger.warning(f"Skipping unsafe path in archive: {member.name}")
            tar.extractall(into, members=safe_members)
    except tarfile.TarError as e:
        raise DataDownloadError(
            f"Downloaded recombination archive from {url} is not a valid tar.gz: {e}"
        ) from e
    return into


def _stage_maps(map_files: list[Path], source: RecombSource, staging: Path) -> None:
    """Write each extracted map into ``staging`` under its canonical name.

    Args:
        map_files: Files the source's ``archive_glob`` matched.
        source: Source whose ``chrom_pattern`` names the chromosome.
        staging: Directory to write into.

    Raises:
        DataDownloadError: If a matched file's name does not carry a
            chromosome, or its content is not a recombination map.
    """
    pattern = re.compile(source.chrom_pattern, re.IGNORECASE)
    for map_file in map_files:
        match = pattern.match(map_file.stem)
        if match is None:
            raise DataDownloadError(
                f"Recombination map {map_file.name!r} does not name a "
                f"chromosome; expected a stem matching "
                f"{source.chrom_pattern!r}. The archive may have changed layout."
            )
        content = ensure_recomb_header(map_file.read_text(), map_file.name)
        (staging / f"chr{match.group(1)}_recomb.tsv").write_text(content)


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
        DataDownloadError: If the download fails, or the archive is corrupt,
            incomplete, or not a recombination map set.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Downloading {source.species} recombination maps...")
    logger.debug(f"Source: {source.url}")

    with tempfile.TemporaryDirectory(dir=output_path.parent) as tmpdir:
        tmp = Path(tmpdir)
        archive = tmp / "maps.tar.gz"
        download_file(source.url, archive, desc="Recombination maps")
        logger.debug(f"Downloaded {archive.stat().st_size / 1024:.1f} KB")

        extracted = _extract_archive(archive, tmp, source.url)
        map_files = sorted(extracted.rglob(source.archive_glob))
        if not map_files:
            logger.error(
                f"Extracted files: {[f.name for f in list(extracted.rglob('*'))[:20]]}"
            )
            raise DataDownloadError(
                f"Could not find chromosome map files matching "
                f"{source.archive_glob!r} in archive from {source.url}"
            )
        logger.debug(f"Found {len(map_files)} chromosome files")

        staging = tmp / "_staging"
        staging.mkdir()
        _stage_maps(map_files, source, staging)
        _publish_map_generation(staging, output_path, source)

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
        output_dir: Directory to save maps. Uses platform cache if None.
        force: Re-download even if files exist.

    Returns:
        Path to the directory containing recombination map files.

    Raises:
        DataDownloadError: If the download fails or the archive is corrupt,
            incomplete, or not a recombination map set.
    """
    output_path = _resolve_map_dir(output_dir)
    if not force and _has_complete_maps(output_path, CANINE_SOURCE):
        return output_path
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
        FileNotFoundError: If map file not found.
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
        raise FileNotFoundError(f"Recombination map not found: {map_file}\n{remedy}")

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
) -> pd.DataFrame:
    """Get recombination rate data for a genomic region.

    Args:
        chrom: Chromosome number.
        start: Start position (bp).
        end: End position (bp).
        species: Species name, alias or record, or None for a caller
            supplying its own maps.
        data_dir: Directory containing recombination maps.
        genome_build: Target genome build (e.g., "canfam4"). If specified and
            different from source data (CanFam3.1), coordinates are lifted over.

    Returns:
        DataFrame with pos and rate columns for the region.

    Note:
        Built-in canine recombination maps are in CanFam3.1 coordinates.
        If genome_build="canfam4", positions are automatically lifted over.
        This requires pyliftover: pip install pyliftover
    """
    record = resolve_species(species)
    df = load_recombination_map(chrom, species=record, data_dir=data_dir)

    source = RECOMB_SOURCES.get(record.key) if record else None
    target_build = assembly_token(genome_build) if genome_build else ""
    if source is not None and target_build in source.liftover_chains:
        logger.debug(f"Lifting over recombination map for chr{chrom} to {target_build}")
        df = liftover_recombination_map(
            df,
            from_build=source.native_build,
            to_build=target_build,
            chrom=chrom,
        )

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
        data_dir: Directory for recombination maps. Uses default if None.

    Returns:
        Path to the recombination maps directory, or None if the species has
        no built-in map set.

    Raises:
        ValidationError: If the species is not one this package knows.
        DataDownloadError: If the species has maps and they could not be
            fetched. Callers that would rather degrade than fail should use
            ``recomb_for_region``, which reports this as a status.
        OSError: If the maps could not be written.
    """
    record = resolve_species(species)
    source = RECOMB_SOURCES.get(record.key) if record else None
    if source is None:
        logger.debug(f"No built-in recombination maps for species: {species}")
        return None

    output_path = _resolve_map_dir(data_dir)

    if _has_complete_maps(output_path, source):
        logger.debug(f"Recombination maps already exist at {output_path}")
        return output_path

    return download_recombination_maps(source, output_path)


class RecombStatus(Enum):
    """Why a region does or does not have recombination rates to draw."""

    OK = "ok"
    NO_MAPS_FOR_SPECIES = "no_maps_for_species"
    NO_MAP_FOR_CHROMOSOME = "no_map_for_chromosome"
    DOWNLOAD_FAILED = "download_failed"
    LIFTOVER_UNAVAILABLE = "liftover_unavailable"


@dataclass(frozen=True)
class RecombResult:
    """The outcome of asking for one region's recombination rates.

    Attributes:
        status: Why there is or is not a frame.
        frame: The region's ``pos`` and ``rate``, set only when status is OK.
        detail: One sentence naming the cause, for the caller to render. Empty
            when status is OK.
    """

    status: RecombStatus
    frame: Optional[pd.DataFrame] = None
    detail: str = ""


def recomb_for_region(
    chrom: int,
    start: int,
    end: int,
    *,
    species: str | Species | None = "canine",
    data_dir: Optional[str] = None,
    genome_build: Optional[str] = None,
) -> RecombResult:
    """Get a region's recombination rates, or say why there are none.

    The one place the "skip the overlay" decision is made. Three layers used
    to make it independently, so whether the user heard about it depended on
    which one fired: a download failure warned, an unsupported species was
    silent, and a missing map file only reached the log. This reports every
    outcome the same way and warns about none of them, leaving the caller to
    render one policy.

    Args:
        chrom: Chromosome number.
        start: Start position (bp).
        end: End position (bp).
        species: Species name, alias or record.
        data_dir: Directory holding the maps. Uses the platform cache if None.
        genome_build: Target build. Coordinates are lifted over if the source
            offers a chain to it.

    Returns:
        A RecombResult carrying the frame, or the status and the reason.

    Raises:
        ValidationError: If the species is not one this package knows.
    """
    record = resolve_species(species)
    try:
        map_dir = ensure_recomb_maps(species=record, data_dir=data_dir)
    except (DataDownloadError, OSError) as e:
        return RecombResult(
            RecombStatus.DOWNLOAD_FAILED,
            detail=f"could not download recombination maps: {e}",
        )

    if map_dir is None:
        name = record.key if record else species
        return RecombResult(
            RecombStatus.NO_MAPS_FOR_SPECIES,
            detail=f"there are no built-in recombination maps for {name!r}",
        )

    try:
        frame = get_recombination_rate_for_region(
            chrom=chrom,
            start=start,
            end=end,
            species=record,
            data_dir=str(map_dir),
            genome_build=genome_build,
        )
    except FileNotFoundError as e:
        return RecombResult(RecombStatus.NO_MAP_FOR_CHROMOSOME, detail=str(e))
    except OptionalDependencyMissing as e:
        return RecombResult(RecombStatus.LIFTOVER_UNAVAILABLE, detail=str(e))

    return RecombResult(RecombStatus.OK, frame=frame)

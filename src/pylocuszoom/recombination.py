"""Recombination rate overlay and data management.

Provides:
- Recombination rate overlay for regional plots
- Download and loading of species-specific recombination maps
- Liftover support for CanFam3.1 to CanFam4 coordinate conversion
"""

import os
import shutil
import tarfile
import tempfile
import uuid
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from ._http import download_file
from ._liftover import PyLiftOverLifter, liftover_positions
from .exceptions import DataDownloadError
from .logging import logger
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

    Attributes:
        species: Species name callers pass, and the key in ``RECOMB_SOURCES``.
        filenames: The complete map set. A directory holding exactly these is
            a cache hit; anything else is downloaded again.
        native_build: Build the published maps are in.
        liftover_chains: Target build token -> chain URL, for the builds the
            maps can be lifted to. A build absent from here is served in
            ``native_build`` coordinates.
        download: Downloads the set, taking ``output_dir`` and ``force``.
    """

    species: str
    filenames: frozenset[str]
    native_build: str
    liftover_chains: dict[str, str]
    download: Callable[..., Path]


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
    except ImportError:
        raise ImportError(
            "pyliftover is required for CanFam4 liftover. "
            "Install it with: pip install pyliftover"
        )

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
    # Determine output directory
    if output_dir is None:
        output_path = get_default_data_dir()
    else:
        output_path = Path(output_dir)

    # Check if already downloaded
    if not force and _has_complete_maps(output_path, CANINE_SOURCE):
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading canine recombination maps from GitHub...")
    logger.debug(f"Source: {CANINE_RECOMB_URL}")

    # Write all per-chromosome files into a staging dir first and only
    # promote them to output_path on full success. A mid-loop raise
    # (e.g. corrupted-archive header rejection) would otherwise leave a
    # partial set of *_recomb.tsv files behind that the cache-hit check
    # at the top of this function might later treat as valid.
    with tempfile.TemporaryDirectory(dir=output_path.parent) as tmpdir:
        # Download tar.gz file with progress bar
        tar_path = Path(tmpdir) / "dog_genetic_maps.tar.gz"

        download_file(
            CANINE_RECOMB_URL,
            tar_path,
            desc="Recombination maps",
        )

        logger.debug(f"Downloaded {tar_path.stat().st_size / 1024:.1f} KB")

        # Extract tar.gz with path traversal protection
        logger.debug("Extracting genetic maps...")
        try:
            with tarfile.open(tar_path, "r:gz") as tar:
                # Filter to prevent path traversal attacks
                safe_members = []
                for member in tar.getmembers():
                    # Resolve the path and ensure it stays within tmpdir
                    member_path = Path(tmpdir) / member.name
                    try:
                        member_path.resolve().relative_to(Path(tmpdir).resolve())
                        safe_members.append(member)
                    except (ValueError, OSError):
                        logger.warning(
                            f"Skipping unsafe path in archive: {member.name}"
                        )
                tar.extractall(tmpdir, members=safe_members)
        except tarfile.TarError as e:
            raise DataDownloadError(
                f"Downloaded recombination archive from {CANINE_RECOMB_URL} "
                f"is not a valid tar.gz: {e}"
            ) from e

        # Find and process the extracted files
        extracted_dir = Path(tmpdir)

        # Look for genetic map files (may be in a subdirectory)
        map_files = list(extracted_dir.rglob("chr*_average_*.txt"))
        if not map_files:
            map_files = list(extracted_dir.rglob("chr*.txt"))
        if not map_files:
            map_files = list(extracted_dir.rglob("*chr*.tsv"))

        if not map_files:
            all_files = list(extracted_dir.rglob("*"))
            logger.error(f"Extracted files: {[f.name for f in all_files[:20]]}")
            raise DataDownloadError("Could not find chromosome map files in archive")

        logger.debug(f"Found {len(map_files)} chromosome files")

        # Stage outputs in tmpdir first; promote to output_path only after
        # every chromosome file is written successfully.
        staging_dir = Path(tmpdir) / "_staging"
        staging_dir.mkdir(exist_ok=True)

        # Copy and rename files
        for map_file in map_files:
            name = map_file.stem
            if "chr" in name.lower():
                chrom = name.lower().split("chr")[-1].split("_")[0].split(".")[0]
                output_file = staging_dir / f"chr{chrom}_recomb.tsv"

                with open(map_file, "r") as f:
                    content = f.read()

                content = ensure_recomb_header(content, map_file.name)

                with open(output_file, "w") as f:
                    f.write(content)

        _publish_map_generation(staging_dir, output_path, CANINE_SOURCE)

    logger.info(f"Recombination maps saved to: {output_path}")
    return output_path


CANINE_SOURCE = RecombSource(
    species="canine",
    filenames=CANINE_MAP_FILENAMES,
    native_build="canfam3",
    liftover_chains={"canfam4": CANFAM3_TO_CANFAM4_CHAIN_URL},
    download=download_canine_recombination_maps,
)

# Species with built-in maps. A species absent from here has no built-in
# source, and its caller supplies data_dir instead.
RECOMB_SOURCES: dict[str, RecombSource] = {CANINE_SOURCE.species: CANINE_SOURCE}


def load_recombination_map(
    chrom: int,
    species: str = "canine",
    data_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Load recombination map for a specific chromosome.

    Args:
        chrom: Chromosome number (1-38 for canine, 1-18 for feline) or 'X'.
        species: Species name ('canine', 'feline').
        data_dir: Directory containing recombination maps.

    Returns:
        DataFrame with columns: pos, rate, cM.

    Raises:
        FileNotFoundError: If map file not found.
    """
    if data_dir is None:
        data_dir = get_default_data_dir()

    data_path = Path(data_dir)
    chrom_str = str(chrom).replace("chr", "")
    map_file = data_path / f"chr{chrom_str}_recomb.tsv"

    if not map_file.exists():
        remedy = (
            f"Run ensure_recomb_maps(species={species!r}) first to download them."
            if species in RECOMB_SOURCES
            else f"There are no built-in recombination maps for {species!r}; "
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
    species: str = "canine",
    data_dir: Optional[str] = None,
    genome_build: Optional[str] = None,
) -> pd.DataFrame:
    """Get recombination rate data for a genomic region.

    Args:
        chrom: Chromosome number.
        start: Start position (bp).
        end: End position (bp).
        species: Species name ('canine', 'feline').
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
    df = load_recombination_map(chrom, species=species, data_dir=data_dir)

    # Liftover if the source publishes a chain to the caller's build
    source = RECOMB_SOURCES.get(species)
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
        chrom_col="",  # Recomb maps don't have chromosome column
        pos_col="pos",
    )

    return region_df[["pos", "rate"]]


def ensure_recomb_maps(
    species: str = "canine",
    data_dir: Optional[str] = None,
) -> Optional[Path]:
    """Ensure recombination maps are available, downloading if needed.

    Args:
        species: Species name ('canine', 'feline', etc.).
        data_dir: Directory for recombination maps. Uses default if None.

    Returns:
        Path to recombination maps directory, or None if species not supported
        or the maps could not be downloaded. A failed download emits a
        UserWarning naming the cause so the missing overlay is not silent.
    """
    source = RECOMB_SOURCES.get(species)
    if source is None:
        logger.debug(f"No built-in recombination maps for species: {species}")
        return None

    if data_dir is not None:
        output_path = Path(data_dir)
    else:
        output_path = get_default_data_dir()

    # Check if maps already exist
    if _has_complete_maps(output_path, source):
        logger.debug(f"Recombination maps already exist at {output_path}")
        return output_path

    # Download maps with error handling
    logger.info(f"Downloading {source.species} recombination maps...")
    try:
        return source.download(output_dir=str(output_path))
    except (DataDownloadError, OSError) as e:
        warnings.warn(
            f"Recombination overlay skipped; could not download recombination "
            f"maps: {e}",
            stacklevel=2,
        )
        return None

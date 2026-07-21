"""Recombination rate overlay and data management.

Provides:
- Recombination rate overlay for regional plots
- Download and loading of species-specific recombination maps
- Liftover support for CanFam3.1 to CanFam4 coordinate conversion
"""

import os
import shutil
import sys
import tarfile
import tempfile
import uuid
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from tqdm import tqdm

from ._liftover import PyLiftOverLifter, liftover_positions
from .logging import logger
from .utils import filter_by_region

# Recombination overlay color
RECOMB_COLOR = "#7FCDFF"  # Light blue
CANINE_MAP_FILENAMES = frozenset(f"chr{chrom}_recomb.tsv" for chrom in range(1, 39))

# Data sources by species
CANINE_RECOMB_URL = (
    "https://github.com/cflerin/dog_recombination/raw/master/dog_genetic_maps.tar.gz"
)

# Liftover chain files
CANFAM3_TO_CANFAM4_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/gbdb/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz"

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
        RuntimeError: If the first token is non-numeric and not a known header.
    """
    lines = content.strip().split("\n")
    first_token = lines[0].split()[0] if lines[0].split() else ""
    normalised_token = first_token.lstrip("#").lower()
    try:
        float(first_token)
        has_header = False
    except ValueError:
        if normalised_token not in _KNOWN_HEADER_TOKENS:
            raise RuntimeError(
                f"Unrecognised first token {first_token!r} in recombination "
                f"map {source_name}; refusing to treat as header. The "
                f"downloaded archive may be corrupted."
            )
        has_header = True
    if not has_header:
        content = _CANONICAL_RECOMB_HEADER + content
    return content


def _normalize_build(build: Optional[str]) -> Optional[str]:
    """Normalize genome build name to canonical form.

    Args:
        build: Build name (e.g., "canfam4", "CanFam4.0", "UU_Cfam_GSD_1.0")

    Returns:
        Normalized build name ("canfam3", "canfam4", or the input lowercased
        for unrecognized builds). None if not specified.
    """
    if build is None:
        return None
    build_lower = build.lower().replace(".", "").replace("_", "")
    if any(x in build_lower for x in ("canfam4", "uucfamgsd")):
        return "canfam4"
    if "canfam3" in build_lower:
        return "canfam3"
    normalized = build.lower()
    logger.debug(
        f"Genome build '{build}' not recognized as canine; using as-is: '{normalized}'"
    )
    return normalized


def get_chain_file_path() -> Path:
    """Get path to the CanFam3 to CanFam4 liftover chain file."""
    return get_default_data_dir() / "canFam3ToCanFam4.over.chain.gz"


def _download_with_progress(
    url: str, dest_path: Path, desc: str = "Downloading"
) -> None:
    """Download a file with a progress bar.

    Args:
        url: URL to download from.
        dest_path: Destination file path.
        desc: Description for the progress bar.
    """
    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))

    with (
        open(dest_path, "wb") as f,
        tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=desc,
            disable=total_size == 0,  # Disable if size unknown
        ) as pbar,
    ):
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))


def download_liftover_chain(force: bool = False) -> Path:
    """Download the CanFam3 to CanFam4 liftover chain file.

    Args:
        force: Re-download even if file exists.

    Returns:
        Path to the downloaded chain file.

    Raises:
        requests.HTTPError: If the download fails with an HTTP error.
    """
    chain_path = get_chain_file_path()

    if chain_path.exists() and not force:
        return chain_path

    chain_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading CanFam3 to CanFam4 liftover chain...")
    logger.debug(f"Source: {CANFAM3_TO_CANFAM4_CHAIN_URL}")

    try:
        _download_with_progress(
            CANFAM3_TO_CANFAM4_CHAIN_URL,
            chain_path,
            desc="Liftover chain",
        )
    except requests.HTTPError as e:
        raise requests.HTTPError(
            f"HTTP error downloading liftover chain file from "
            f"{CANFAM3_TO_CANFAM4_CHAIN_URL}: {e}"
        ) from e

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

    Returns platform-appropriate cache directory:
    - macOS/Linux: ~/.cache/pylocuszoom (or $XDG_CACHE_HOME if set)
    - Windows: %LOCALAPPDATA%/pylocuszoom
    - Databricks: /dbfs/FileStore/reference_data/recombination_maps
    """
    if sys.platform == "win32":  # Windows
        base = Path(os.environ.get("LOCALAPPDATA", Path.home()))
    elif os.path.exists("/dbfs"):  # Databricks
        return Path("/dbfs/FileStore/reference_data/recombination_maps")
    else:
        # macOS and Linux
        xdg_cache = os.environ.get("XDG_CACHE_HOME")
        if xdg_cache:
            base = Path(xdg_cache)
        else:
            base = Path.home() / ".cache"

    return base / "pylocuszoom" / "recombination_maps"


def _has_complete_canine_maps(path: Path) -> bool:
    """Return whether path contains exactly the expected canine map set."""
    if not path.exists():
        return False
    present = {map_path.name for map_path in path.glob("chr*_recomb.tsv")}
    return present == CANINE_MAP_FILENAMES


def _publish_map_generation(staging_dir: Path, output_path: Path) -> Path:
    """Publish a complete map generation behind an atomically replaced link."""
    parent = output_path.parent
    parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    generation = parent / f".{output_path.name}.generation-{token}"

    if output_path.exists():
        shutil.copytree(output_path, generation)
        for old_map in generation.glob("chr*_recomb.tsv"):
            old_map.unlink()
    else:
        generation.mkdir()

    for staged_file in staging_dir.glob("chr*_recomb.tsv"):
        shutil.move(str(staged_file), generation / staged_file.name)

    if not _has_complete_canine_maps(generation):
        shutil.rmtree(generation)
        raise RuntimeError(
            "Downloaded recombination archive does not contain the complete canine "
            "map set (chromosomes 1-38)"
        )

    pending_link = parent / f".{output_path.name}.pending-{token}"
    pending_link.symlink_to(generation.name, target_is_directory=True)

    if output_path.is_symlink() or not output_path.exists():
        previous_generation = (
            output_path.resolve() if output_path.is_symlink() else None
        )
        os.replace(pending_link, output_path)
    else:
        previous_generation = None
        legacy = parent / f".{output_path.name}.legacy-{token}"
        os.replace(output_path, legacy)
        try:
            os.replace(pending_link, output_path)
        except BaseException:
            os.replace(legacy, output_path)
            shutil.rmtree(generation, ignore_errors=True)
            pending_link.unlink(missing_ok=True)
            raise
        shutil.rmtree(legacy)

    if (
        previous_generation is not None
        and previous_generation != generation
        and previous_generation.parent == parent
        and previous_generation.name.startswith(f".{output_path.name}.generation-")
    ):
        shutil.rmtree(previous_generation, ignore_errors=True)

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
    """
    # Determine output directory
    if output_dir is None:
        output_path = get_default_data_dir()
    else:
        output_path = Path(output_dir)

    # Check if already downloaded
    if not force and _has_complete_canine_maps(output_path):
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

        _download_with_progress(
            CANINE_RECOMB_URL,
            tar_path,
            desc="Recombination maps",
        )

        logger.debug(f"Downloaded {tar_path.stat().st_size / 1024:.1f} KB")

        # Extract tar.gz with path traversal protection
        logger.debug("Extracting genetic maps...")
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
                    logger.warning(f"Skipping unsafe path in archive: {member.name}")
            tar.extractall(tmpdir, members=safe_members)

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
            raise RuntimeError("Could not find chromosome map files in archive")

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

        _publish_map_generation(staging_dir, output_path)

    logger.info(f"Recombination maps saved to: {output_path}")
    return output_path


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
        raise FileNotFoundError(
            f"Recombination map not found: {map_file}\n"
            f"Run download_{species}_recombination_maps() first to download the data."
        )

    df = pd.read_csv(map_file, sep="\t")

    for col in ["pos", "rate"]:
        original = df[col]
        df[col] = pd.to_numeric(df[col], errors="coerce")
        bad_mask = df[col].isna() & original.notna()
        bad_count = bad_mask.sum()
        if bad_count > 0:
            sample_vals = original[bad_mask].head(3).tolist()
            logger.warning(
                f"Recombination map chr{chrom_str}: {bad_count} non-numeric "
                f"values in '{col}' column dropped. Sample values: {sample_vals}"
            )
    if "cM" in df.columns:
        original = df["cM"]
        df["cM"] = pd.to_numeric(df["cM"], errors="coerce")
        bad_mask = df["cM"].isna() & original.notna()
        bad_count = bad_mask.sum()
        if bad_count > 0:
            sample_vals = original[bad_mask].head(3).tolist()
            logger.warning(
                f"Recombination map chr{chrom_str}: {bad_count} non-numeric "
                f"values in 'cM' column dropped. Sample values: {sample_vals}"
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

    # Liftover if needed
    build = _normalize_build(genome_build)
    if species == "canine" and build == "canfam4":
        logger.debug(f"Lifting over recombination map for chr{chrom} to CanFam4")
        df = liftover_recombination_map(
            df, from_build="canfam3", to_build="canfam4", chrom=chrom
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
        or download fails.
    """
    if species != "canine":
        logger.debug(f"No built-in recombination maps for species: {species}")
        return None

    if data_dir is not None:
        output_path = Path(data_dir)
    else:
        output_path = get_default_data_dir()

    # Check if maps already exist
    if _has_complete_canine_maps(output_path):
        logger.debug(f"Recombination maps already exist at {output_path}")
        return output_path

    # Download maps with error handling
    logger.info("Downloading canine recombination maps...")
    try:
        return download_canine_recombination_maps(output_dir=str(output_path))
    except (requests.RequestException, tarfile.TarError, OSError) as e:
        logger.error(f"Could not download recombination maps: {e!r}")
        return None

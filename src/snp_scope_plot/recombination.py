"""Recombination rate overlay and data management.

Provides:
- Recombination rate overlay for regional plots
- Download and loading of species-specific recombination maps
- Liftover support for CanFam3.1 to CanFam4 coordinate conversion
"""

import os
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Optional

import pandas as pd
from matplotlib.axes import Axes

from .logging import logger

# Recombination overlay color
RECOMB_COLOR = "#7FCDFF"  # Light blue

# Data sources by species
DOG_RECOMB_URL = (
    "https://github.com/cflerin/dog_recombination/raw/master/dog_genetic_maps.tar.gz"
)

# Liftover chain files
CANFAM3_TO_CANFAM4_CHAIN_URL = "https://hgdownload.soe.ucsc.edu/gbdb/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz"


def _normalize_build(build: Optional[str]) -> Optional[str]:
    """Normalize genome build name to canonical form.

    Args:
        build: Build name (e.g., "canfam4", "CanFam4.0", "UU_Cfam_GSD_1.0")

    Returns:
        Normalized build name ("canfam3" or "canfam4"), or None if not specified.
    """
    if build is None:
        return None
    build_lower = build.lower().replace(".", "").replace("_", "")
    if "canfam4" in build_lower or "uucfamgsd" in build_lower:
        return "canfam4"
    if "canfam3" in build_lower:
        return "canfam3"
    return build.lower()


def get_chain_file_path() -> Path:
    """Get path to the CanFam3 to CanFam4 liftover chain file."""
    return get_default_data_dir() / "canFam3ToCanFam4.over.chain.gz"


def download_liftover_chain(force: bool = False) -> Path:
    """Download the CanFam3 to CanFam4 liftover chain file.

    Args:
        force: Re-download even if file exists.

    Returns:
        Path to the downloaded chain file.
    """
    chain_path = get_chain_file_path()

    if chain_path.exists() and not force:
        return chain_path

    chain_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading CanFam3 to CanFam4 liftover chain...")
    logger.debug(f"Source: {CANFAM3_TO_CANFAM4_CHAIN_URL}")

    try:
        urllib.request.urlretrieve(CANFAM3_TO_CANFAM4_CHAIN_URL, chain_path)
    except Exception as e:
        logger.debug(f"urllib download failed: {e}")
        try:
            import requests

            response = requests.get(CANFAM3_TO_CANFAM4_CHAIN_URL, timeout=60)
            response.raise_for_status()
            chain_path.write_bytes(response.content)
        except ImportError:
            raise RuntimeError(
                "Failed to download. Install requests: pip install requests"
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
    from pyliftover import LiftOver

    # Download chain file if needed
    chain_path = download_liftover_chain()

    logger.debug(f"Lifting over coordinates from {from_build} to {to_build}")
    lo = LiftOver(str(chain_path))

    # Get chromosome for each position
    if "chr" in recomb_df.columns:
        chroms = recomb_df["chr"].astype(str)
    elif chrom is not None:
        chroms = pd.Series([str(chrom)] * len(recomb_df))
    else:
        raise ValueError("Either 'chr' column or chrom parameter required")

    # Liftover each position
    new_positions = []
    keep_mask = []

    for chr_val, pos in zip(chroms, recomb_df["pos"]):
        chr_str = f"chr{chr_val}" if not str(chr_val).startswith("chr") else chr_val
        result = lo.convert_coordinate(chr_str, int(pos))

        if result and len(result) > 0:
            # Take first mapping (usually the only one)
            _, new_pos, _, _ = result[0]
            new_positions.append(int(new_pos))
            keep_mask.append(True)
        else:
            new_positions.append(None)
            keep_mask.append(False)

    # Create output DataFrame
    result_df = recomb_df.copy()
    result_df["pos"] = new_positions
    result_df = result_df[keep_mask].copy()

    unmapped = len(recomb_df) - len(result_df)
    if unmapped > 0:
        logger.debug(f"Dropped {unmapped} positions that failed to liftover")

    return result_df.sort_values("pos").reset_index(drop=True)


def get_default_data_dir() -> Path:
    """Get default directory for recombination map data.

    Returns platform-appropriate cache directory:
    - macOS: ~/Library/Caches/snp-scope-plot
    - Linux: ~/.cache/snp-scope-plot
    - Windows: %LOCALAPPDATA%/snp-scope-plot
    """
    if os.name == "nt":  # Windows
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

    return base / "snp-scope-plot" / "recombination_maps"


def download_dog_recombination_maps(
    output_dir: Optional[str] = None,
    force: bool = False,
) -> Path:
    """Download dog recombination rate maps from Campbell et al. 2016.

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
    if output_path.exists() and not force:
        existing_files = list(output_path.glob("chr*_recomb.tsv"))
        if len(existing_files) >= 38:  # 38 autosomes + X
            return output_path

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading dog recombination maps from GitHub...")
    logger.debug(f"Source: {DOG_RECOMB_URL}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Download tar.gz file
        tar_path = Path(tmpdir) / "dog_genetic_maps.tar.gz"

        try:
            urllib.request.urlretrieve(DOG_RECOMB_URL, tar_path)
        except Exception as e:
            logger.debug(f"urllib download failed: {e}")
            logger.debug("Trying alternative method with requests...")
            try:
                import requests

                response = requests.get(DOG_RECOMB_URL, timeout=60)
                response.raise_for_status()
                tar_path.write_bytes(response.content)
            except ImportError:
                raise RuntimeError(
                    "Failed to download. Install requests: pip install requests"
                )

        logger.debug(f"Downloaded {tar_path.stat().st_size / 1024:.1f} KB")

        # Extract tar.gz
        logger.debug("Extracting genetic maps...")
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(tmpdir)

        # Find and process the extracted files
        extracted_dir = Path(tmpdir)

        # Look for genetic map files (may be in a subdirectory)
        map_files = list(extracted_dir.rglob("chr*.txt"))
        if not map_files:
            map_files = list(extracted_dir.rglob("*chr*.tsv"))

        if not map_files:
            all_files = list(extracted_dir.rglob("*"))
            logger.error(f"Extracted files: {[f.name for f in all_files[:20]]}")
            raise RuntimeError("Could not find chromosome map files in archive")

        logger.debug(f"Found {len(map_files)} chromosome files")

        # Copy and rename files
        for map_file in map_files:
            name = map_file.stem
            if "chr" in name.lower():
                chrom = name.lower().split("chr")[-1].split("_")[0].split(".")[0]
                output_file = output_path / f"chr{chrom}_recomb.tsv"

                with open(map_file, "r") as f:
                    content = f.read()

                # Ensure header is present
                lines = content.strip().split("\n")
                if not lines[0].startswith("chr") and not lines[0].startswith("pos"):
                    content = "chr\tpos\trate\tcM\n" + content

                with open(output_file, "w") as f:
                    f.write(content)

    logger.info(f"Recombination maps saved to: {output_path}")
    return output_path


def load_recombination_map(
    chrom: int,
    species: str = "dog",
    data_dir: Optional[str] = None,
) -> pd.DataFrame:
    """Load recombination map for a specific chromosome.

    Args:
        chrom: Chromosome number (1-38 for dog, 1-18 for cat) or 'X'.
        species: Species name ('dog', 'cat').
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

    # Ensure numeric columns
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["rate"] = pd.to_numeric(df["rate"], errors="coerce")
    if "cM" in df.columns:
        df["cM"] = pd.to_numeric(df["cM"], errors="coerce")

    return df.dropna(subset=["pos", "rate"])


def get_recombination_rate_for_region(
    chrom: int,
    start: int,
    end: int,
    species: str = "dog",
    data_dir: Optional[str] = None,
    genome_build: Optional[str] = None,
) -> pd.DataFrame:
    """Get recombination rate data for a genomic region.

    Args:
        chrom: Chromosome number.
        start: Start position (bp).
        end: End position (bp).
        species: Species name ('dog', 'cat').
        data_dir: Directory containing recombination maps.
        genome_build: Target genome build (e.g., "canfam4"). If specified and
            different from source data (CanFam3.1), coordinates are lifted over.

    Returns:
        DataFrame with pos and rate columns for the region.

    Note:
        Built-in dog recombination maps are in CanFam3.1 coordinates.
        If genome_build="canfam4", positions are automatically lifted over.
        This requires pyliftover: pip install pyliftover
    """
    df = load_recombination_map(chrom, species=species, data_dir=data_dir)

    # Liftover if needed
    build = _normalize_build(genome_build)
    if species == "dog" and build == "canfam4":
        logger.debug(f"Lifting over recombination map for chr{chrom} to CanFam4")
        df = liftover_recombination_map(
            df, from_build="canfam3", to_build="canfam4", chrom=chrom
        )

    # Filter to region
    region_df = df[(df["pos"] >= start) & (df["pos"] <= end)].copy()

    return region_df[["pos", "rate"]]


def add_recombination_overlay(
    ax: Axes,
    recomb_df: pd.DataFrame,
    start: int,
    end: int,
) -> Optional[Axes]:
    """Add recombination rate as secondary y-axis overlay.

    Plots recombination rate (cM/Mb) as a light blue line on a
    secondary y-axis, styled to match LocusZoom.

    Args:
        ax: Primary matplotlib axes object.
        recomb_df: DataFrame with 'pos' and 'rate' columns.
        start: Region start position.
        end: Region end position.

    Returns:
        Secondary axes object for recombination rate, or None if no data.
    """
    # Create secondary y-axis
    recomb_ax = ax.twinx()

    # Filter to region
    region_recomb = recomb_df[
        (recomb_df["pos"] >= start) & (recomb_df["pos"] <= end)
    ].copy()

    if region_recomb.empty:
        recomb_ax.set_visible(False)
        return None

    # Plot recombination rate as light blue line
    recomb_ax.plot(
        region_recomb["pos"],
        region_recomb["rate"],
        color=RECOMB_COLOR,
        linewidth=1.5,
        alpha=0.7,
        zorder=0,  # Behind scatter points
    )

    # Fill under curve
    recomb_ax.fill_between(
        region_recomb["pos"],
        0,
        region_recomb["rate"],
        color=RECOMB_COLOR,
        alpha=0.15,
        zorder=0,
    )

    # Format secondary axis
    recomb_ax.set_ylabel("Recombination rate (cM/Mb)", color=RECOMB_COLOR, fontsize=9)
    recomb_ax.tick_params(axis="y", labelcolor=RECOMB_COLOR, labelsize=8)
    recomb_ax.set_ylim(bottom=0)

    # Don't let recomb rate overwhelm the plot
    max_rate = region_recomb["rate"].max()
    recomb_ax.set_ylim(0, max(max_rate * 1.2, 20))

    # Remove top spine for cleaner look
    recomb_ax.spines["top"].set_visible(False)

    return recomb_ax

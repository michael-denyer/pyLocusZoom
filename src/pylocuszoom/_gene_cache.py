# src/pylocuszoom/_gene_cache.py
"""Disk cache for fetched gene annotations, shared by the reference sources.

Both ``ensembl.py`` and ``ucsc.py`` cache the same shape of result under the
same rules, so the key derivation and the CSV round-trip live here rather than
once per client. Each source owns its own cache root, so a region fetched from
Ensembl and the same region fetched from UCSC never collide.
"""

import hashlib
from pathlib import Path

import pandas as pd

from .exceptions import ValidationError
from .logging import logger
from .utils import _platform_cache_base, normalize_chrom


def cache_root(source: str) -> Path:
    """Get the cache directory for one reference source.

    Shares the platform cache base with recombination maps (see
    ``utils._platform_cache_base``), so ``$XDG_CACHE_HOME`` is honored on macOS
    and Linux and Databricks routes to ``/dbfs/FileStore/reference_data``.

    Args:
        source: Source leaf directory, e.g. ``"ensembl"`` or ``"ucsc"``.

    Returns:
        Path to the cache directory (created if it doesn't exist).
    """
    path = _platform_cache_base() / source
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_species_dir(cache_dir: Path, species: str) -> Path:
    """Resolve a species subdirectory and validate it stays within cache_dir.

    Prevents path traversal from untrusted species strings (e.g. ``"../../etc"``
    would escape the cache root).

    Args:
        cache_dir: Root cache directory.
        species: Resolved species name, already mapped through its source's
            alias table.

    Returns:
        Resolved Path to the species subdirectory.

    Raises:
        ValidationError: If the resolved path escapes cache_dir.
    """
    species_dir = (cache_dir / species).resolve()
    if not species_dir.is_relative_to(cache_dir.resolve()):
        raise ValidationError(
            f"Invalid species name: {species!r} (resolved path escapes cache directory)"
        )
    return species_dir


def cache_key(
    species: str,
    chrom: str,
    start: int,
    end: int,
    build_token: str = "",
) -> str:
    """Generate the cache key for a region.

    The build token is part of the key so two plots of the same region under
    different builds never share an entry. Including it also orphans every
    entry written before builds were tracked, whose assembly is unknowable.
    """
    key_str = f"{species}_{chrom}_{start}_{end}_{build_token}"
    return hashlib.md5(key_str.encode()).hexdigest()[:16]


def load_genes(
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    build_token: str = "",
) -> pd.DataFrame | None:
    """Load a cached gene frame, or None on a miss or an unreadable file."""
    key = cache_key(species, normalize_chrom(chrom), start, end, build_token)
    cache_file = safe_species_dir(cache_dir, species) / f"genes_{key}.csv"

    if not cache_file.exists():
        return None

    try:
        logger.debug(f"Cache hit: {cache_file}")
        return pd.read_csv(cache_file)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        # EmptyDataError subclasses ValueError, not ParserError, so it needs
        # naming explicitly. Releases before 2.1.1 wrote column-less CSVs.
        logger.warning(f"Corrupt cache file {cache_file}, ignoring: {e}")
        return None


def save_genes(
    df: pd.DataFrame,
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    build_token: str = "",
) -> None:
    """Write a gene frame to the cache, logging rather than raising on failure."""
    key = cache_key(species, normalize_chrom(chrom), start, end, build_token)
    species_dir = safe_species_dir(cache_dir, species)
    species_dir.mkdir(parents=True, exist_ok=True)

    cache_file = species_dir / f"genes_{key}.csv"
    try:
        df.to_csv(cache_file, index=False)
        logger.debug(f"Cached genes to: {cache_file}")
    except OSError as e:
        logger.warning(f"Failed to write gene cache {cache_file}: {e}")


def clear_cache(cache_dir: Path, species: str | None = None) -> int:
    """Delete cached gene files, optionally for one species only.

    Args:
        cache_dir: Cache directory to clear.
        species: If given, only this resolved species subdirectory is cleared.

    Returns:
        Number of files deleted.
    """
    if not cache_dir.exists():
        return 0

    if species is not None:
        search_dirs = [safe_species_dir(cache_dir, species)]
    else:
        search_dirs = [d for d in cache_dir.iterdir() if d.is_dir()]

    deleted = 0
    for directory in search_dirs:
        if not directory.exists():
            continue
        for cache_file in directory.glob("*.csv"):
            try:
                cache_file.unlink()
                deleted += 1
            except OSError as e:
                logger.warning(f"Failed to delete {cache_file}: {e}")
    return deleted

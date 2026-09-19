# src/pylocuszoom/_gene_cache.py
"""Disk cache for fetched gene annotations, shared by the reference sources.

Both ``ensembl.py`` and ``ucsc.py`` cache the same shape of result under the
same rules, so key derivation and the CSV archive round-trip live here
rather than once per client. Each source owns its own cache root, so a region fetched from
Ensembl and the same region fetched from UCSC never collide.

One archive contains both frames and is published with one replacement. Legacy
CSV pairs are misses because they cannot prove both frames belong together.
"""

import hashlib
import os
import tempfile
from pathlib import Path
from zipfile import BadZipFile, ZipFile
from zlib import error as ZlibError

import pandas as pd

from ._gene_source import GeneAnnotations
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


def _entry_file(
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    build_token: str,
) -> Path:
    """Resolve the single archive holding a complete cache entry."""
    key = cache_key(species, normalize_chrom(chrom), start, end, build_token)
    species_dir = safe_species_dir(cache_dir, species)
    return species_dir / f"annotations_{key}.zip"


def load_annotations(
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    build_token: str = "",
) -> GeneAnnotations | None:
    """Load both frames from one published archive, or return None on a miss."""
    entry = _entry_file(cache_dir, species, chrom, start, end, build_token)
    try:
        with ZipFile(entry) as archive:
            with archive.open("genes.csv") as genes, archive.open("exons.csv") as exons:
                result = GeneAnnotations(pd.read_csv(genes), pd.read_csv(exons))
        logger.debug(f"Cache hit: {entry}")
        return result
    except FileNotFoundError:
        return None
    except (
        OSError,
        BadZipFile,
        EOFError,
        RuntimeError,  # Encrypted entries and unsupported ZIP compression.
        ZlibError,
        KeyError,
        UnicodeDecodeError,
        pd.errors.ParserError,
        pd.errors.EmptyDataError,
    ) as e:
        logger.warning(f"Corrupt cache file for {entry}, ignoring: {e}")
        return None


def save_annotations(
    annotations: GeneAnnotations,
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    build_token: str = "",
) -> None:
    """Publish a complete entry, leaving the old entry intact on failure."""
    entry = _entry_file(cache_dir, species, chrom, start, end, build_token)
    partial_path = None
    try:
        entry.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=entry.parent, prefix=f".{entry.stem}.", suffix=".part", delete=False
        ) as partial:
            partial_path = Path(partial.name)
        with ZipFile(partial_path, "w") as archive:
            with archive.open("genes.csv", "w") as genes:
                annotations.genes.to_csv(genes, index=False)
            with archive.open("exons.csv", "w") as exons:
                annotations.exons.to_csv(exons, index=False)
        os.replace(partial_path, entry)
        logger.debug(f"Cached annotations to: {entry}")
    except OSError as e:
        logger.warning(f"Failed to write gene cache {entry}: {e}")
    finally:
        if partial_path is not None:
            try:
                partial_path.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(
                    f"Failed to clean up gene cache staging {partial_path}: {e}"
                )


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
        for cache_file in (
            *directory.glob("*.csv"),
            *directory.glob("annotations_*.zip"),
        ):
            try:
                cache_file.unlink()
                deleted += 1
            except OSError as e:
                logger.warning(f"Failed to delete {cache_file}: {e}")
    return deleted

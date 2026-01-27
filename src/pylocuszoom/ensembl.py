# src/pylocuszoom/ensembl.py
"""Ensembl REST API integration for reference data fetching.

Provides functions to fetch gene and exon annotations from the Ensembl REST API
(https://rest.ensembl.org) for any species.

Note: Recombination rates are NOT available from Ensembl for most species.
Use species-specific recombination maps instead (see recombination.py).
"""

import hashlib
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests

from .logging import logger
from .utils import ValidationError

# Ensembl API limits regions to 5Mb
ENSEMBL_MAX_REGION_SIZE = 5_000_000

# Species name aliases -> Ensembl species names
SPECIES_ALIASES: dict[str, str] = {
    # Canine
    "canine": "canis_lupus_familiaris",
    "dog": "canis_lupus_familiaris",
    "canis_familiaris": "canis_lupus_familiaris",
    # Feline
    "feline": "felis_catus",
    "cat": "felis_catus",
    # Human
    "human": "homo_sapiens",
    # Mouse
    "mouse": "mus_musculus",
    # Rat
    "rat": "rattus_norvegicus",
}


ENSEMBL_REST_URL = "https://rest.ensembl.org"
ENSEMBL_REQUEST_TIMEOUT = 30  # seconds
ENSEMBL_MAX_RETRIES = 3
ENSEMBL_RETRY_DELAY = 1.0  # seconds, doubles on each retry


def _normalize_chrom(chrom: str | int) -> str:
    """Normalize chromosome name by removing 'chr' prefix."""
    return str(chrom).replace("chr", "")


def _validate_region_size(start: int, end: int, context: str) -> None:
    """Validate region size is within Ensembl API limits.

    Args:
        start: Region start position.
        end: Region end position.
        context: Context for error message (e.g., "genes_df", "exons_df").

    Raises:
        ValidationError: If region exceeds 5Mb limit.
    """
    region_size = end - start
    if region_size > ENSEMBL_MAX_REGION_SIZE:
        raise ValidationError(
            f"Region size {region_size:,} bp exceeds Ensembl API limit of 5Mb. "
            f"Please use a smaller region or provide {context} directly."
        )


def get_ensembl_species_name(species: str) -> str:
    """Convert species alias to Ensembl species name.

    Args:
        species: Species name or alias (e.g., "canine", "dog", "human").

    Returns:
        Ensembl-compatible species name (e.g., "canis_lupus_familiaris").
    """
    return SPECIES_ALIASES.get(species.lower(), species.lower())


def get_ensembl_cache_dir() -> Path:
    """Get the cache directory for Ensembl data.

    Uses same base location as recombination maps: ~/.cache/snp-scope-plot/ensembl

    Returns:
        Path to cache directory (created if doesn't exist).
    """
    if sys.platform == "darwin":
        base = Path.home() / ".cache"
    elif sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))

    cache_dir = base / "snp-scope-plot" / "ensembl"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _cache_key(species: str, chrom: str, start: int, end: int) -> str:
    """Generate cache key for a region."""
    key_str = f"{species}_{chrom}_{start}_{end}"
    return hashlib.md5(key_str.encode()).hexdigest()[:16]


def get_cached_genes(
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
) -> pd.DataFrame | None:
    """Load cached genes if available.

    Args:
        cache_dir: Cache directory path.
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position.
        end: Region end position.

    Returns:
        DataFrame if cache hit, None if cache miss.
    """
    ensembl_species = get_ensembl_species_name(species)
    chrom_str = _normalize_chrom(chrom)
    cache_key = _cache_key(ensembl_species, chrom_str, start, end)

    species_dir = cache_dir / ensembl_species
    cache_file = species_dir / f"genes_{cache_key}.csv"

    if not cache_file.exists():
        return None

    logger.debug(f"Cache hit: {cache_file}")
    return pd.read_csv(cache_file)


def save_cached_genes(
    df: pd.DataFrame,
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
) -> None:
    """Save genes to cache as CSV.

    Args:
        df: DataFrame with gene annotations to cache.
        cache_dir: Cache directory path.
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position.
        end: Region end position.
    """
    ensembl_species = get_ensembl_species_name(species)
    chrom_str = _normalize_chrom(chrom)
    cache_key = _cache_key(ensembl_species, chrom_str, start, end)

    species_dir = cache_dir / ensembl_species
    species_dir.mkdir(parents=True, exist_ok=True)

    cache_file = species_dir / f"genes_{cache_key}.csv"
    df.to_csv(cache_file, index=False)
    logger.debug(f"Cached genes to: {cache_file}")


def _make_ensembl_request(
    url: str,
    params: dict,
    max_retries: int = ENSEMBL_MAX_RETRIES,
    raise_on_error: bool = False,
) -> list | None:
    """Make request to Ensembl API with retry logic.

    Args:
        url: API endpoint URL.
        params: Query parameters.
        max_retries: Maximum retry attempts for retryable errors.
        raise_on_error: If True, raise exception on error instead of returning None.

    Returns:
        JSON response as list, or None on non-retryable error.

    Raises:
        ValidationError: If raise_on_error=True and request fails.
    """
    delay = ENSEMBL_RETRY_DELAY

    for attempt in range(max_retries):
        try:
            response = requests.get(
                url,
                params=params,
                headers={"Content-Type": "application/json"},
                timeout=ENSEMBL_REQUEST_TIMEOUT,
            )
        except requests.RequestException as e:
            logger.warning(f"Ensembl API request failed (attempt {attempt + 1}): {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
                delay *= 2
                continue
            if raise_on_error:
                raise ValidationError(
                    f"Ensembl API request failed after {max_retries} attempts: {e}"
                )
            return None

        # Success
        if response.ok:
            return response.json()

        # Retryable errors (429 rate limit, 503 service unavailable)
        if response.status_code in (429, 503) and attempt < max_retries - 1:
            logger.warning(
                f"Ensembl API returned {response.status_code} "
                f"(attempt {attempt + 1}), retrying..."
            )
            time.sleep(delay)
            delay *= 2
            continue

        # Non-retryable error
        error_msg = f"Ensembl API error {response.status_code}: {response.text[:200]}"
        logger.warning(error_msg)
        if raise_on_error:
            raise ValidationError(error_msg)
        return None

    return None


def fetch_genes_from_ensembl(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    biotype: str = "protein_coding",
    raise_on_error: bool = False,
) -> pd.DataFrame:
    """Fetch gene annotations from Ensembl REST API.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        biotype: Gene biotype filter (default: protein_coding).
        raise_on_error: If True, raise ValidationError on API errors.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand, gene_id, biotype.
        Returns empty DataFrame on API error (unless raise_on_error=True).

    Raises:
        ValidationError: If region > 5Mb or if raise_on_error=True and API fails.
    """
    _validate_region_size(start, end, "genes_df")

    ensembl_species = get_ensembl_species_name(species)
    chrom_str = _normalize_chrom(chrom)

    # Build region string
    region = f"{chrom_str}:{start}-{end}"

    # Build API URL
    url = f"{ENSEMBL_REST_URL}/overlap/region/{ensembl_species}/{region}"
    params = {"feature": "gene", "biotype": biotype}

    logger.debug(f"Fetching genes from Ensembl: {url}")

    data = _make_ensembl_request(url, params, raise_on_error=raise_on_error)

    if data is None:
        return pd.DataFrame()

    if not data:
        logger.debug(f"No genes found in region {region}")
        return pd.DataFrame()

    # Convert to DataFrame
    records = []
    for gene in data:
        if gene.get("feature_type") != "gene":
            continue
        records.append(
            {
                "chr": str(gene.get("seq_region_name", chrom_str)),
                "start": gene.get("start"),
                "end": gene.get("end"),
                "gene_name": gene.get("external_name", gene.get("id", "")),
                "strand": "+" if gene.get("strand", 1) == 1 else "-",
                "gene_id": gene.get("id", ""),
                "biotype": gene.get("biotype", ""),
            }
        )

    df = pd.DataFrame(records)
    logger.debug(f"Fetched {len(df)} genes from Ensembl")
    return df


def fetch_exons_from_ensembl(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    raise_on_error: bool = False,
) -> pd.DataFrame:
    """Fetch exon annotations from Ensembl REST API.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        raise_on_error: If True, raise ValidationError on API errors.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, exon_id, transcript_id.
        Returns empty DataFrame on API error (unless raise_on_error=True).

    Raises:
        ValidationError: If region > 5Mb or if raise_on_error=True and API fails.
    """
    _validate_region_size(start, end, "exons_df")

    ensembl_species = get_ensembl_species_name(species)
    chrom_str = _normalize_chrom(chrom)
    region = f"{chrom_str}:{start}-{end}"

    url = f"{ENSEMBL_REST_URL}/overlap/region/{ensembl_species}/{region}"
    params = {"feature": "exon"}

    logger.debug(f"Fetching exons from Ensembl: {url}")

    data = _make_ensembl_request(url, params, raise_on_error=raise_on_error)

    if data is None:
        return pd.DataFrame()

    if not data:
        return pd.DataFrame()

    records = []
    for exon in data:
        if exon.get("feature_type") != "exon":
            continue
        records.append(
            {
                "chr": str(exon.get("seq_region_name", chrom_str)),
                "start": exon.get("start"),
                "end": exon.get("end"),
                "gene_name": "",  # Exon endpoint doesn't include gene name
                "exon_id": exon.get("id", ""),
                "transcript_id": exon.get("Parent", ""),
            }
        )

    df = pd.DataFrame(records)
    logger.debug(f"Fetched {len(df)} exons from Ensembl")
    return df


def get_genes_for_region(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    cache_dir: Path | None = None,
    use_cache: bool = True,
    include_exons: bool = False,
    raise_on_error: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Get gene annotations for a genomic region.

    Checks cache first, fetches from Ensembl API if not cached.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        cache_dir: Cache directory (uses default if None).
        use_cache: Whether to use disk cache.
        include_exons: If True, also fetch exons and return tuple (genes_df, exons_df).
        raise_on_error: If True, raise ValidationError on API errors.

    Returns:
        If include_exons=False: DataFrame with gene annotations.
        If include_exons=True: Tuple of (genes_df, exons_df).

    Raises:
        ValidationError: If region > 5Mb or if raise_on_error=True and API fails.

    Note:
        Gene annotations are cached to disk. Exons are fetched from the API
        on each call when include_exons=True (not cached separately).
    """
    if cache_dir is None:
        cache_dir = get_ensembl_cache_dir()

    chrom_str = _normalize_chrom(chrom)

    # Check cache first
    if use_cache:
        cached = get_cached_genes(cache_dir, species, chrom_str, start, end)
        if cached is not None:
            if include_exons:
                # Exons not cached separately (yet)
                exons_df = fetch_exons_from_ensembl(
                    species, chrom_str, start, end, raise_on_error=raise_on_error
                )
                return cached, exons_df
            return cached

    # Fetch from Ensembl API
    genes_df = fetch_genes_from_ensembl(
        species, chrom_str, start, end, raise_on_error=raise_on_error
    )

    # Cache the result (even if empty, to avoid repeated API calls for gene-sparse regions)
    if use_cache:
        save_cached_genes(genes_df, cache_dir, species, chrom_str, start, end)

    if include_exons:
        exons_df = fetch_exons_from_ensembl(
            species, chrom_str, start, end, raise_on_error=raise_on_error
        )
        return genes_df, exons_df

    return genes_df


def clear_ensembl_cache(
    cache_dir: Path | None = None,
    species: str | None = None,
) -> int:
    """Clear cached Ensembl data.

    Args:
        cache_dir: Cache directory (uses default if None).
        species: If provided, only clear cache for this species.

    Returns:
        Number of files deleted.
    """
    if cache_dir is None:
        cache_dir = get_ensembl_cache_dir()

    deleted = 0

    if species:
        # Clear only specific species
        ensembl_species = get_ensembl_species_name(species)
        species_dir = cache_dir / ensembl_species
        if species_dir.exists():
            for cache_file in species_dir.glob("*.csv"):
                cache_file.unlink()
                deleted += 1
    else:
        # Clear all species
        for cache_file in cache_dir.glob("**/*.csv"):
            cache_file.unlink()
            deleted += 1

    logger.info(f"Cleared {deleted} cached Ensembl files from {cache_dir}")
    return deleted

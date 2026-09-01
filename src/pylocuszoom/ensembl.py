# src/pylocuszoom/ensembl.py
"""Ensembl REST API integration for reference data fetching.

Provides functions to fetch gene and exon annotations from the Ensembl REST API
(https://rest.ensembl.org) for any species.

Ensembl serves exactly one reference assembly per species and answers a request
naming any other with that same assembly and an HTTP 200, so the caller's
genome build cannot be honoured and a mismatch is invisible without checking
``assembly_name`` on the response. Every fetch here therefore takes an optional
``genome_build`` and warns when the two disagree. Release 116 was the last on
this REST platform, so retired assemblies such as CanFam3.1 and FelCat9 will
not reappear on it.

Note: Recombination rates are NOT available from Ensembl for most species.
Use species-specific recombination maps instead (see recombination.py).
"""

import warnings
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from ._gene_cache import cache_root, clear_cache, load_genes, save_genes
from ._http import request_json
from .exceptions import EnsemblAPIError, ValidationError
from .logging import logger
from .utils import assembly_token, normalize_chrom

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


def _response_assembly(features: list) -> str:
    """Read the assembly Ensembl actually served from an overlap response."""
    for feature in features:
        assembly = feature.get("assembly_name")
        if assembly:
            return str(assembly)
    return ""


def _warn_on_assembly_mismatch(
    assembly: str, genome_build: str | None, species: str
) -> None:
    """Warn when Ensembl's assembly differs from the caller's genome build.

    Ensembl serves exactly one reference assembly per species and silently
    ignores a ``coord_system_version`` asking for any other, so a mismatch
    yields plausible-looking genes in the wrong coordinate system rather than
    an error. Dog is the worst case in practice: Ensembl retired CanFam3.1 and
    now serves ROS_Cfam_1.0 only, which puts ATP9B at chr1:938,796 instead of
    chr1:1,136,865, a shift of roughly 198 kb.
    """
    if not assembly or not genome_build:
        return
    if assembly_token(assembly) == assembly_token(genome_build):
        return
    warnings.warn(
        f"Ensembl returned {species} annotations on assembly {assembly!r}, but "
        f"genome_build is {genome_build!r}. Ensembl serves one assembly per "
        f"species and ignores requests for any other, so these gene coordinates "
        f"do not line up with your data. Supply genes_df yourself in "
        f"{genome_build!r} coordinates, or set genome_build={assembly!r}.",
        UserWarning,
        stacklevel=2,
    )


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

    Returns:
        Path to cache directory (created if doesn't exist).
    """
    return cache_root("ensembl")


def get_cached_genes(
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    genome_build: str | None = None,
) -> pd.DataFrame | None:
    """Load cached genes if available.

    Args:
        cache_dir: Cache directory path.
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position.
        end: Region end position.
        genome_build: Build the caller's data is in; part of the cache key and
            checked against the assembly recorded in the cached frame.

    Returns:
        DataFrame if cache hit, None if cache miss.
    """
    ensembl_species = get_ensembl_species_name(species)
    df = load_genes(
        cache_dir,
        ensembl_species,
        chrom,
        start,
        end,
        build_token=assembly_token(genome_build or ""),
    )
    if df is None:
        return None

    if "assembly" in df.columns and not df.empty:
        _warn_on_assembly_mismatch(
            str(df["assembly"].iloc[0]), genome_build, ensembl_species
        )
    return df


def save_cached_genes(
    df: pd.DataFrame,
    cache_dir: Path,
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    genome_build: str | None = None,
) -> None:
    """Save genes to cache as CSV.

    Args:
        df: DataFrame with gene annotations to cache.
        cache_dir: Cache directory path.
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position.
        end: Region end position.
        genome_build: Build the caller's data is in; part of the cache key.
    """
    save_genes(
        df,
        cache_dir,
        get_ensembl_species_name(species),
        chrom,
        start,
        end,
        build_token=assembly_token(genome_build or ""),
    )


def _gene_record(feature: dict, chrom_str: str) -> dict:
    """Build a gene DataFrame row from an Ensembl overlap feature."""
    return {
        "chr": str(feature.get("seq_region_name", chrom_str)),
        "start": feature.get("start"),
        "end": feature.get("end"),
        "gene_name": feature.get("external_name", feature.get("id", "")),
        "strand": "+" if feature.get("strand", 1) == 1 else "-",
        "gene_id": feature.get("id", ""),
        "biotype": feature.get("biotype", ""),
        "assembly": str(feature.get("assembly_name", "")),
    }


def _exon_record(feature: dict, chrom_str: str) -> dict:
    """Build an exon DataFrame row from an Ensembl overlap feature."""
    return {
        "chr": str(feature.get("seq_region_name", chrom_str)),
        "start": feature.get("start"),
        "end": feature.get("end"),
        "gene_name": "",  # Exon endpoint doesn't include gene name
        "exon_id": feature.get("id", ""),
        "transcript_id": feature.get("Parent", ""),
        "assembly": str(feature.get("assembly_name", "")),
    }


def _empty_feature_frame(
    record_builder: Callable[[dict, str], dict], chrom_str: str
) -> pd.DataFrame:
    """Build an empty frame carrying the record builder's columns.

    A bare ``pd.DataFrame()`` has no columns, so it serialises to a one-byte CSV
    that ``pd.read_csv`` cannot parse back. Deriving the columns from the record
    builder's own keys keeps the schema single-sourced.
    """
    return pd.DataFrame(columns=list(record_builder({}, chrom_str)))


def _fetch_overlap_features(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    *,
    params: dict,
    feature_type: str,
    record_builder: Callable[[dict, str], dict],
    raise_on_error: bool = False,
    genome_build: str | None = None,
) -> pd.DataFrame:
    """Fetch features from the Ensembl overlap/region endpoint.

    Shared by fetch_genes_from_ensembl and fetch_exons_from_ensembl, which differ
    only in the request params, the feature_type they keep, and how each feature
    maps to a DataFrame row (record_builder).

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        params: Query parameters for the overlap request.
        feature_type: Ensembl feature_type to keep (e.g. "gene", "exon").
        record_builder: Maps a kept feature and chrom to a DataFrame row.
        raise_on_error: If True, raise EnsemblAPIError on API errors.
        genome_build: Build the caller's data is in; warns if Ensembl serves a
            different assembly.

    Returns:
        DataFrame of built records, always carrying the record_builder's
        columns; empty on API error or an empty region.

    Raises:
        ValidationError: If region > 5Mb.
        EnsemblAPIError: If raise_on_error=True and the API fails.
    """
    _validate_region_size(start, end, f"{feature_type}s_df")

    ensembl_species = get_ensembl_species_name(species)
    chrom_str = normalize_chrom(chrom)
    region = f"{chrom_str}:{start}-{end}"
    url = f"{ENSEMBL_REST_URL}/overlap/region/{ensembl_species}/{region}"

    logger.debug(f"Fetching {feature_type}s from Ensembl: {url}")

    empty = _empty_feature_frame(record_builder, chrom_str)

    try:
        data = request_json(
            url,
            params,
            error_cls=EnsemblAPIError,
            service="Ensembl",
            headers={"Content-Type": "application/json"},
            timeout=ENSEMBL_REQUEST_TIMEOUT,
            max_retries=ENSEMBL_MAX_RETRIES,
            retry_delay=ENSEMBL_RETRY_DELAY,
        )
    except EnsemblAPIError:
        if raise_on_error:
            raise
        return empty

    if not data:
        logger.debug(f"No {feature_type}s found in region {region}")
        return empty

    _warn_on_assembly_mismatch(_response_assembly(data), genome_build, ensembl_species)

    records = [
        record_builder(feature, chrom_str)
        for feature in data
        if feature.get("feature_type") == feature_type
    ]
    if not records:
        logger.debug(f"No {feature_type}s found in region {region}")
        return empty

    df = pd.DataFrame(records)
    logger.debug(f"Fetched {len(df)} {feature_type}s from Ensembl")
    return df


def fetch_genes_from_ensembl(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    biotype: str = "protein_coding",
    raise_on_error: bool = False,
    genome_build: str | None = None,
) -> pd.DataFrame:
    """Fetch gene annotations from Ensembl REST API.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        biotype: Gene biotype filter (default: protein_coding).
        raise_on_error: If True, raise ValidationError on API errors.
        genome_build: Build the caller's data is in; warns if Ensembl serves a
            different assembly.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand, gene_id,
        biotype, assembly.
        Returns empty DataFrame on API error (unless raise_on_error=True).

    Raises:
        ValidationError: If region > 5Mb or if raise_on_error=True and API fails.
    """
    return _fetch_overlap_features(
        species,
        chrom,
        start,
        end,
        params={"feature": "gene", "biotype": biotype},
        feature_type="gene",
        record_builder=_gene_record,
        raise_on_error=raise_on_error,
        genome_build=genome_build,
    )


def fetch_exons_from_ensembl(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    raise_on_error: bool = False,
    genome_build: str | None = None,
) -> pd.DataFrame:
    """Fetch exon annotations from Ensembl REST API.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        raise_on_error: If True, raise ValidationError on API errors.
        genome_build: Build the caller's data is in; warns if Ensembl serves a
            different assembly.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, exon_id,
        transcript_id, assembly.
        Returns empty DataFrame on API error (unless raise_on_error=True).

    Raises:
        ValidationError: If region > 5Mb or if raise_on_error=True and API fails.
    """
    return _fetch_overlap_features(
        species,
        chrom,
        start,
        end,
        params={"feature": "exon"},
        feature_type="exon",
        record_builder=_exon_record,
        raise_on_error=raise_on_error,
        genome_build=genome_build,
    )


def get_genes_for_region(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    cache_dir: Path | None = None,
    use_cache: bool = True,
    include_exons: bool = False,
    raise_on_error: bool = False,
    genome_build: str | None = None,
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
        genome_build: Build the caller's data is in. Part of the cache key, and
            warns when Ensembl serves annotations on a different assembly.

    Returns:
        If include_exons=False: DataFrame with gene annotations.
        If include_exons=True: Tuple of (genes_df, exons_df).
        Both carry an ``assembly`` column naming the assembly Ensembl served.

    Raises:
        ValidationError: If region > 5Mb or if raise_on_error=True and API fails.

    Note:
        Gene annotations are cached to disk. Exons are fetched from the API
        on each call when include_exons=True (not cached separately).
    """
    if cache_dir is None:
        cache_dir = get_ensembl_cache_dir()

    chrom_str = normalize_chrom(chrom)

    # Check cache first
    if use_cache:
        cached = get_cached_genes(
            cache_dir, species, chrom_str, start, end, genome_build=genome_build
        )
        if cached is not None:
            if include_exons:
                # Exons not cached separately (yet)
                exons_df = fetch_exons_from_ensembl(
                    species,
                    chrom_str,
                    start,
                    end,
                    raise_on_error=raise_on_error,
                    genome_build=genome_build,
                )
                return cached, exons_df
            return cached

    # Fetch from Ensembl API. Ask it to raise so a service failure stays
    # distinguishable from a region that genuinely has no genes; only the
    # latter is safe to cache.
    fetch_failed = False
    try:
        genes_df = fetch_genes_from_ensembl(
            species,
            chrom_str,
            start,
            end,
            raise_on_error=True,
            genome_build=genome_build,
        )
    except EnsemblAPIError:
        if raise_on_error:
            raise
        fetch_failed = True
        genes_df = _empty_feature_frame(_gene_record, chrom_str)

    # Cache the result (even if empty, to avoid repeated API calls for gene-sparse
    # regions). A failed fetch is never cached: it would be indistinguishable from
    # an empty region on reload and would permanently hide the region's genes.
    if use_cache and not fetch_failed:
        save_cached_genes(
            genes_df,
            cache_dir,
            species,
            chrom_str,
            start,
            end,
            genome_build=genome_build,
        )

    if include_exons:
        exons_df = fetch_exons_from_ensembl(
            species,
            chrom_str,
            start,
            end,
            raise_on_error=raise_on_error,
            genome_build=genome_build,
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

    ensembl_species = get_ensembl_species_name(species) if species else None
    deleted = clear_cache(cache_dir, ensembl_species)
    logger.info(f"Cleared {deleted} cached Ensembl files from {cache_dir}")
    return deleted

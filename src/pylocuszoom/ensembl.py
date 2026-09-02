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
from pathlib import Path

import pandas as pd

from ._gene_cache import cache_root
from ._http import request_json
from .exceptions import EnsemblAPIError, ValidationError
from .logging import logger
from .reference_genes import (
    EXON_COLUMNS,
    GENE_COLUMNS,
    clear_gene_cache,
    ensembl_source,
    get_genes_for_build,
)
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


def warn_on_cached_assembly(
    cached: pd.DataFrame, genome_build: str | None, species: str
) -> None:
    """Repeat the assembly-mismatch warning for a frame reloaded from cache.

    The warning fires on the original fetch, but a later session reads the same
    genes off disk without asking Ensembl, so without this the mismatch would
    be announced exactly once in the cache entry's lifetime.
    """
    if cached.empty or "assembly" not in cached.columns:
        return
    _warn_on_assembly_mismatch(str(cached["assembly"].iloc[0]), genome_build, species)


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


def _exon_record(feature: dict, chrom_str: str, gene_name: str) -> dict:
    """Build an exon DataFrame row from an Ensembl overlap feature."""
    return {
        "chr": str(feature.get("seq_region_name", chrom_str)),
        "start": feature.get("start"),
        "end": feature.get("end"),
        "gene_name": gene_name,
        "exon_id": feature.get("id", ""),
        "transcript_id": feature.get("Parent", ""),
        "assembly": str(feature.get("assembly_name", "")),
    }


def _frame(records: list[dict], columns: tuple[str, ...]) -> pd.DataFrame:
    """Build a frame that carries ``columns`` even when there are no records."""
    return pd.DataFrame(records) if records else pd.DataFrame(columns=list(columns))


def _by_feature_type(data: list[dict]) -> dict[str, list[dict]]:
    """Split an overlap response into its gene, transcript and exon features."""
    buckets: dict[str, list[dict]] = {"gene": [], "transcript": [], "exon": []}
    for feature in data:
        bucket = buckets.get(str(feature.get("feature_type")))
        if bucket is not None:
            bucket.append(feature)
    return buckets


def fetch_overlap_frames(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    genome_build: str | None = None,
    biotype: str = "protein_coding",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch the gene and exon frames for a region in one overlap request.

    Exon features name only their transcript, and transcript features name
    only their gene, so asking for all three feature types in one request is
    what lets an exon carry the gene symbol the gene track matches it by. An
    exon whose gene is absent from the response (its biotype was filtered out)
    is dropped rather than returned unattached.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        genome_build: Build the caller's data is in; warns if Ensembl serves a
            different assembly.
        biotype: Gene biotype filter.

    Returns:
        Tuple of (genes_df, exons_df), each carrying its full column set even
        when the region has nothing of that kind.

    Raises:
        ValidationError: If region > 5Mb.
        EnsemblAPIError: If the API fails.
    """
    _validate_region_size(start, end, "genes_df")

    ensembl_species = get_ensembl_species_name(species)
    chrom_str = normalize_chrom(chrom)
    region = f"{chrom_str}:{start}-{end}"
    url = f"{ENSEMBL_REST_URL}/overlap/region/{ensembl_species}/{region}"

    logger.debug(f"Fetching genes and exons from Ensembl: {url}")

    data = request_json(
        url,
        {"feature": ["gene", "transcript", "exon"], "biotype": biotype},
        error_cls=EnsemblAPIError,
        service="Ensembl",
        headers={"Content-Type": "application/json"},
        timeout=ENSEMBL_REQUEST_TIMEOUT,
        max_retries=ENSEMBL_MAX_RETRIES,
        retry_delay=ENSEMBL_RETRY_DELAY,
    )

    if not data:
        logger.debug(f"No genes found in region {region}")
        return _frame([], GENE_COLUMNS), _frame([], EXON_COLUMNS)

    _warn_on_assembly_mismatch(_response_assembly(data), genome_build, ensembl_species)

    features = _by_feature_type(data)
    genes = [_gene_record(feature, chrom_str) for feature in features["gene"]]
    symbol_of_gene = {
        feature.get("id"): record["gene_name"]
        for feature, record in zip(features["gene"], genes)
    }
    gene_of_transcript = {
        feature.get("id"): feature.get("Parent") for feature in features["transcript"]
    }
    exons = []
    for feature in features["exon"]:
        symbol = symbol_of_gene.get(gene_of_transcript.get(feature.get("Parent")), "")
        if symbol:
            exons.append(_exon_record(feature, chrom_str, symbol))

    logger.debug(f"Fetched {len(genes)} genes and {len(exons)} exons from Ensembl")
    return _frame(genes, GENE_COLUMNS), _frame(exons, EXON_COLUMNS)


def _frames_or_empty(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    raise_on_error: bool,
    genome_build: str | None = None,
    biotype: str = "protein_coding",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch frames, translating a service failure to empties if allowed."""
    try:
        return fetch_overlap_frames(species, chrom, start, end, genome_build, biotype)
    except EnsemblAPIError:
        if raise_on_error:
            raise
        return _frame([], GENE_COLUMNS), _frame([], EXON_COLUMNS)


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
        raise_on_error: If True, raise EnsemblAPIError on API errors.
        genome_build: Build the caller's data is in; warns if Ensembl serves a
            different assembly.

    Returns:
        DataFrame with the columns in ``reference_genes.GENE_COLUMNS``.
        Returns empty DataFrame on API error (unless raise_on_error=True).

    Raises:
        ValidationError: If region > 5Mb.
        EnsemblAPIError: If raise_on_error=True and the API fails.
    """
    return _frames_or_empty(
        species, chrom, start, end, raise_on_error, genome_build, biotype
    )[0]


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
        raise_on_error: If True, raise EnsemblAPIError on API errors.
        genome_build: Build the caller's data is in; warns if Ensembl serves a
            different assembly.

    Returns:
        DataFrame with the columns in ``reference_genes.EXON_COLUMNS``.
        Returns empty DataFrame on API error (unless raise_on_error=True).

    Raises:
        ValidationError: If region > 5Mb.
        EnsemblAPIError: If raise_on_error=True and the API fails.
    """
    return _frames_or_empty(species, chrom, start, end, raise_on_error, genome_build)[1]


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
    """Get gene annotations for a genomic region from Ensembl.

    Checks cache first, fetches from the Ensembl API if not cached. Bypasses
    the build-to-source routing in ``reference_genes``, so it always asks
    Ensembl even for a build only UCSC can serve.

    Args:
        species: Species name or alias.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        cache_dir: Cache directory (uses default if None).
        use_cache: Whether to use disk cache.
        include_exons: If True, also fetch exons and return tuple (genes_df, exons_df).
        raise_on_error: If True, raise EnsemblAPIError on API errors.
        genome_build: Build the caller's data is in. Part of the cache key, and
            warns when Ensembl serves annotations on a different assembly.

    Returns:
        If include_exons=False: DataFrame with gene annotations.
        If include_exons=True: Tuple of (genes_df, exons_df).
        Both carry an ``assembly`` column naming the assembly Ensembl served.

    Raises:
        ValidationError: If region > 5Mb.
        EnsemblAPIError: If raise_on_error=True and the API fails.

    Note:
        Gene annotations are cached to disk. Exons are fetched from the API
        on each call when include_exons=True (not cached separately).
    """
    return get_genes_for_build(
        species,
        chrom,
        start,
        end,
        genome_build=genome_build,
        cache_dir=cache_dir,
        use_cache=use_cache,
        include_exons=include_exons,
        raise_on_error=raise_on_error,
        source=ensembl_source(species, genome_build),
    )


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
    return clear_gene_cache(
        "ensembl",
        cache_dir,
        get_ensembl_species_name(species) if species else None,
    )

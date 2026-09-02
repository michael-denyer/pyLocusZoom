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

import pandas as pd

from ._gene_source import (
    EXON_COLUMNS,
    GENE_COLUMNS,
    GeneAnnotations,
    GeneSource,
    empty_annotations,
    empty_frame,
)
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


def _validate_region_size(start: int, end: int) -> None:
    """Validate region size is within Ensembl API limits.

    Args:
        start: Region start position.
        end: Region end position.

    Raises:
        ValidationError: If region exceeds 5Mb limit.
    """
    region_size = end - start
    if region_size > ENSEMBL_MAX_REGION_SIZE:
        raise ValidationError(
            f"Region size {region_size:,} bp exceeds Ensembl API limit of 5Mb. "
            f"Please use a smaller region or provide genes_df directly."
        )


def get_ensembl_species_name(species: str) -> str:
    """Convert species alias to Ensembl species name.

    Args:
        species: Species name or alias (e.g., "canine", "dog", "human").

    Returns:
        Ensembl-compatible species name (e.g., "canis_lupus_familiaris").
    """
    return SPECIES_ALIASES.get(species.lower(), species.lower())


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
) -> GeneAnnotations:
    """Fetch the genes and exons for a region in one overlap request.

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
        The region's annotations.

    Raises:
        ValidationError: If region > 5Mb.
        EnsemblAPIError: If the API fails.
    """
    _validate_region_size(start, end)

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
        return empty_annotations()

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
    return GeneAnnotations(
        pd.DataFrame(genes) if genes else empty_frame(GENE_COLUMNS),
        pd.DataFrame(exons) if exons else empty_frame(EXON_COLUMNS),
    )


def ensembl_source(species: str, genome_build: str | None = None) -> GeneSource:
    """Build the GeneSource for one species on Ensembl's current assembly."""
    ensembl_species = get_ensembl_species_name(species)
    return GeneSource(
        name="ensembl",
        cache_species=ensembl_species,
        build_token=assembly_token(genome_build or ""),
        fetch=lambda chrom, start, end: fetch_overlap_frames(
            species, chrom, start, end, genome_build=genome_build
        ),
        on_cache_hit=lambda cached: warn_on_cached_assembly(
            cached, genome_build, ensembl_species
        ),
    )

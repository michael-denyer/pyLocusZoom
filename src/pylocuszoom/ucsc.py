# src/pylocuszoom/ucsc.py
"""UCSC REST API client for gene annotations on assemblies Ensembl has retired.

Ensembl serves one reference assembly per species and release 116 was the last
on the legacy REST platform, so CanFam3.1, CanFam4 and FelCat9 have no Ensembl
source at any URL; the archive REST hosts redirect to a help page rather than
serving data. UCSC still hosts all three, which is why gene tracks for those
builds come from here.

The track is ``ncbiRefSeq`` rather than ``ensGene`` because UCSC's ensGene
carries only stable IDs (ENSCAFG..., ENSFCAG...) in place of gene symbols, and
a gene track labelled with accession numbers is unreadable.

Frames match the shape ``ensembl.py`` returns, including the ``assembly``
column, so either source drops into the same plot.
"""

from pathlib import Path

import pandas as pd

from ._gene_cache import cache_root, clear_cache, load_genes, save_genes
from ._http import request_json
from .exceptions import UCSCAPIError
from .logging import logger
from .utils import normalize_chrom

UCSC_REST_URL = "https://api.genome.ucsc.edu"
UCSC_REQUEST_TIMEOUT = 30  # seconds
UCSC_MAX_RETRIES = 3
UCSC_RETRY_DELAY = 1.0  # seconds, doubles on each retry
UCSC_GENE_TRACK = "ncbiRefSeq"

# UCSC genome name -> the assembly name recorded on every returned row.
# Every genome reachable through reference_genes.UCSC_BUILDS belongs here.
UCSC_ASSEMBLY_NAMES: dict[str, str] = {
    "canFam3": "CanFam3.1",
    "canFam4": "UU_Cfam_GSD_1.0",
    "felCat9": "Felis_catus_9.0",
}

# RefSeq accession prefixes for transcripts that code for protein.
_CODING_PREFIXES = ("NM_", "XM_")


def get_ucsc_cache_dir() -> Path:
    """Get the cache directory for UCSC gene data."""
    return cache_root("ucsc")


def _fetch_track(
    ucsc_genome: str, chrom: str | int, start: int, end: int
) -> list[dict]:
    """Fetch raw ncbiRefSeq transcript rows overlapping a region."""
    chrom_str = normalize_chrom(chrom)
    payload = request_json(
        f"{UCSC_REST_URL}/getData/track",
        {
            "genome": ucsc_genome,
            "track": UCSC_GENE_TRACK,
            "chrom": f"chr{chrom_str}",
            "start": start,
            "end": end,
        },
        error_cls=UCSCAPIError,
        service="UCSC",
        timeout=UCSC_REQUEST_TIMEOUT,
        max_retries=UCSC_MAX_RETRIES,
        retry_delay=UCSC_RETRY_DELAY,
    )
    rows = payload.get(UCSC_GENE_TRACK, [])
    # UCSC returns a chrom-keyed dict when the request spans the whole genome.
    if isinstance(rows, dict):
        rows = rows.get(f"chr{chrom_str}", [])
    return rows


def _is_coding(row: dict) -> bool:
    """A gene codes for protein when any transcript accession is NM_/XM_."""
    return str(row.get("name", "")).startswith(_CODING_PREFIXES)


def _gene_record(row: dict, chrom_str: str, assembly: str) -> dict:
    """Seed a gene row from the first transcript carrying its symbol."""
    symbol = row.get("name2") or row.get("name", "")
    return {
        "chr": chrom_str,
        "start": int(row.get("txStart", 0)) + 1,
        "end": int(row.get("txEnd", 0)),
        "gene_name": symbol,
        "strand": row.get("strand", "+"),
        "gene_id": symbol,
        "biotype": "protein_coding" if _is_coding(row) else "non_coding",
        "assembly": assembly,
    }


def _exon_record(
    row: dict,
    chrom_str: str,
    assembly: str,
    exon_start: str = "0",
    exon_end: str = "0",
    index: int = 0,
) -> dict:
    """Build one exon row from a transcript row and one exon's coordinates."""
    transcript = str(row.get("name", ""))
    return {
        "chr": chrom_str,
        "start": int(exon_start) + 1,
        "end": int(exon_end),
        "gene_name": row.get("name2", ""),
        "exon_id": f"{transcript}_exon{index + 1}",
        "transcript_id": transcript,
        "assembly": assembly,
    }


def _empty_frame(record_builder) -> pd.DataFrame:
    """Build an empty frame carrying the record builder's columns.

    A bare ``pd.DataFrame()`` has no columns, so it serialises to a one-byte
    CSV that ``pd.read_csv`` cannot parse back. Deriving the columns from the
    builder keeps the schema single-sourced.
    """
    return pd.DataFrame(columns=list(record_builder({}, "", "")))


def fetch_genes_from_ucsc(
    ucsc_genome: str,
    chrom: str | int,
    start: int,
    end: int,
    biotype: str = "protein_coding",
    raise_on_error: bool = False,
) -> pd.DataFrame:
    """Fetch gene annotations from UCSC, collapsed from transcripts to genes.

    ncbiRefSeq is a transcript-level track, so the many transcripts sharing a
    symbol are collapsed into one row spanning the widest of them. UCSC
    coordinates are 0-based half-open and are converted to the 1-based
    inclusive convention Ensembl and the rest of pyLocusZoom use.

    Args:
        ucsc_genome: UCSC genome name, e.g. ``"canFam3"``.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        biotype: Gene biotype filter. ncbiRefSeq carries no finer biotypes
            than the accession prefix, so the only meaningful values are
            ``"protein_coding"`` (at least one NM_/XM_ transcript),
            ``"non_coding"`` (none), and None or "" to keep everything; any
            other value matches nothing.
        raise_on_error: If True, raise UCSCAPIError on API errors.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand, gene_id,
        biotype, assembly. Empty on API error unless raise_on_error=True.

    Raises:
        UCSCAPIError: If raise_on_error=True and the API fails.
    """
    empty = _empty_frame(_gene_record)
    try:
        rows = _fetch_track(ucsc_genome, chrom, start, end)
    except UCSCAPIError:
        if raise_on_error:
            raise
        return empty

    assembly = UCSC_ASSEMBLY_NAMES.get(ucsc_genome, ucsc_genome)
    chrom_str = normalize_chrom(chrom)

    genes: dict[str, dict] = {}
    for row in rows:
        symbol = row.get("name2") or row.get("name")
        if not symbol:
            continue
        gene = genes.get(symbol)
        if gene is None:
            genes[symbol] = _gene_record(row, chrom_str, assembly)
            continue
        gene["start"] = min(gene["start"], int(row["txStart"]) + 1)
        gene["end"] = max(gene["end"], int(row["txEnd"]))
        if _is_coding(row):
            gene["biotype"] = "protein_coding"

    records = list(genes.values())
    if biotype:
        records = [g for g in records if g["biotype"] == biotype]

    if not records:
        logger.debug(f"No genes found in {ucsc_genome} {chrom_str}:{start}-{end}")
        return empty

    logger.debug(f"Fetched {len(records)} genes from UCSC {ucsc_genome}")
    return pd.DataFrame(records)


def fetch_exons_from_ucsc(
    ucsc_genome: str,
    chrom: str | int,
    start: int,
    end: int,
    raise_on_error: bool = False,
) -> pd.DataFrame:
    """Fetch exon annotations from UCSC.

    Every transcript's exons are returned, matching what Ensembl's overlap
    endpoint gives for the same region rather than picking one transcript.

    Args:
        ucsc_genome: UCSC genome name, e.g. ``"canFam3"``.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        raise_on_error: If True, raise UCSCAPIError on API errors.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, exon_id,
        transcript_id, assembly. Empty on API error unless raise_on_error=True.

    Raises:
        UCSCAPIError: If raise_on_error=True and the API fails.
    """
    empty = _empty_frame(_exon_record)
    try:
        rows = _fetch_track(ucsc_genome, chrom, start, end)
    except UCSCAPIError:
        if raise_on_error:
            raise
        return empty

    assembly = UCSC_ASSEMBLY_NAMES.get(ucsc_genome, ucsc_genome)
    chrom_str = normalize_chrom(chrom)

    records = []
    for row in rows:
        starts = str(row.get("exonStarts", "")).strip(",")
        ends = str(row.get("exonEnds", "")).strip(",")
        if not starts or not ends:
            continue
        records.extend(
            _exon_record(row, chrom_str, assembly, exon_start, exon_end, index)
            for index, (exon_start, exon_end) in enumerate(
                zip(starts.split(","), ends.split(","))
            )
        )

    if not records:
        logger.debug(f"No exons found in {ucsc_genome} {chrom_str}:{start}-{end}")
        return empty

    logger.debug(f"Fetched {len(records)} exons from UCSC {ucsc_genome}")
    return pd.DataFrame(records)


def get_genes_for_region_ucsc(
    ucsc_genome: str,
    chrom: str | int,
    start: int,
    end: int,
    cache_dir: Path | None = None,
    use_cache: bool = True,
    include_exons: bool = False,
    raise_on_error: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Get gene annotations for a region from UCSC, with disk caching.

    Mirrors ``ensembl.get_genes_for_region``: the cache is keyed by genome and
    region, an empty region is cached so gene-sparse regions stop re-requesting,
    and a failed fetch is never cached because it would be indistinguishable
    from an empty region on reload.

    Args:
        ucsc_genome: UCSC genome name, e.g. ``"canFam3"``.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        cache_dir: Cache directory (uses default if None).
        use_cache: Whether to use the disk cache.
        include_exons: If True, return a (genes_df, exons_df) tuple.
        raise_on_error: If True, raise UCSCAPIError on API errors.

    Returns:
        Gene annotations, or a (genes_df, exons_df) tuple when include_exons.

    Raises:
        UCSCAPIError: If raise_on_error=True and the API fails.
    """
    if cache_dir is None:
        cache_dir = get_ucsc_cache_dir()

    chrom_str = normalize_chrom(chrom)

    def exons():
        return fetch_exons_from_ucsc(
            ucsc_genome, chrom_str, start, end, raise_on_error=raise_on_error
        )

    if use_cache:
        cached = load_genes(cache_dir, ucsc_genome, chrom_str, start, end)
        if cached is not None:
            return (cached, exons()) if include_exons else cached

    fetch_failed = False
    try:
        genes_df = fetch_genes_from_ucsc(
            ucsc_genome, chrom_str, start, end, raise_on_error=True
        )
    except UCSCAPIError:
        if raise_on_error:
            raise
        fetch_failed = True
        genes_df = _empty_frame(_gene_record)

    if use_cache and not fetch_failed:
        save_genes(genes_df, cache_dir, ucsc_genome, chrom_str, start, end)

    return (genes_df, exons()) if include_exons else genes_df


def clear_ucsc_cache(
    cache_dir: Path | None = None, ucsc_genome: str | None = None
) -> int:
    """Clear cached UCSC gene data.

    Args:
        cache_dir: Cache directory (uses default if None).
        ucsc_genome: If given, only clear this genome's cache.

    Returns:
        Number of files deleted.
    """
    if cache_dir is None:
        cache_dir = get_ucsc_cache_dir()
    deleted = clear_cache(cache_dir, ucsc_genome)
    logger.info(f"Cleared {deleted} cached UCSC files from {cache_dir}")
    return deleted

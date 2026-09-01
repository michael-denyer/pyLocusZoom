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

import time

import pandas as pd
import requests

from ._gene_cache import cache_root, clear_cache, load_genes, save_genes
from .exceptions import UCSCAPIError
from .logging import logger
from .utils import normalize_chrom

UCSC_REST_URL = "https://api.genome.ucsc.edu"
UCSC_REQUEST_TIMEOUT = 30  # seconds
UCSC_MAX_RETRIES = 3
UCSC_RETRY_DELAY = 1.0  # seconds, doubles on each retry
UCSC_GENE_TRACK = "ncbiRefSeq"

# UCSC genome name -> the assembly name recorded on every returned row.
UCSC_ASSEMBLY_NAMES: dict[str, str] = {
    "canFam3": "CanFam3.1",
    "canFam4": "UU_Cfam_GSD_1.0",
    "canFam6": "Dog10K_Boxer_Tasha",
    "felCat9": "Felis_catus_9.0",
}

# RefSeq accession prefixes for transcripts that code for protein.
_CODING_PREFIXES = ("NM_", "XM_")


def get_ucsc_cache_dir():
    """Get the cache directory for UCSC gene data."""
    return cache_root("ucsc")


def _make_ucsc_request(url: str, params: dict) -> dict:
    """Make a request to the UCSC API with retry on rate limits and outages.

    Always raises on failure; callers that want an empty result instead
    translate UCSCAPIError at their boundary.

    Raises:
        UCSCAPIError: If the request ultimately fails.
    """
    delay = UCSC_RETRY_DELAY

    for attempt in range(UCSC_MAX_RETRIES):
        try:
            response = requests.get(url, params=params, timeout=UCSC_REQUEST_TIMEOUT)
        except requests.RequestException as e:
            logger.warning(f"UCSC API request failed (attempt {attempt + 1}): {e}")
            if attempt < UCSC_MAX_RETRIES - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise UCSCAPIError(
                f"UCSC API request failed after {UCSC_MAX_RETRIES} attempts: {e}"
            )

        if response.ok:
            try:
                return response.json()
            except (ValueError, requests.exceptions.JSONDecodeError) as e:
                raise UCSCAPIError(f"UCSC API returned invalid JSON: {e}")

        if response.status_code in (429, 503) and attempt < UCSC_MAX_RETRIES - 1:
            logger.warning(
                f"UCSC API returned {response.status_code} "
                f"(attempt {attempt + 1}), retrying..."
            )
            time.sleep(delay)
            delay *= 2
            continue

        raise UCSCAPIError(
            f"UCSC API error {response.status_code}: {response.text[:200]}"
        )

    raise UCSCAPIError(
        f"UCSC API request failed after {UCSC_MAX_RETRIES} attempts (rate limited)"
    )


def _fetch_track(
    ucsc_genome: str, chrom: str | int, start: int, end: int
) -> list[dict]:
    """Fetch raw ncbiRefSeq transcript rows overlapping a region."""
    chrom_str = normalize_chrom(chrom)
    payload = _make_ucsc_request(
        f"{UCSC_REST_URL}/getData/track",
        {
            "genome": ucsc_genome,
            "track": UCSC_GENE_TRACK,
            "chrom": f"chr{chrom_str}",
            "start": start,
            "end": end,
        },
    )
    rows = payload.get(UCSC_GENE_TRACK, [])
    # UCSC returns a chrom-keyed dict when the request spans the whole genome.
    if isinstance(rows, dict):
        rows = rows.get(f"chr{chrom_str}", [])
    return rows


def _assembly_for(ucsc_genome: str) -> str:
    return UCSC_ASSEMBLY_NAMES.get(ucsc_genome, ucsc_genome)


def _gene_columns() -> list[str]:
    return [
        "chr",
        "start",
        "end",
        "gene_name",
        "strand",
        "gene_id",
        "biotype",
        "assembly",
    ]


def _exon_columns() -> list[str]:
    return ["chr", "start", "end", "gene_name", "exon_id", "transcript_id", "assembly"]


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
        biotype: Gene biotype filter. ``"protein_coding"`` keeps genes with at
            least one NM_/XM_ transcript; None or "" keeps everything.
        raise_on_error: If True, raise UCSCAPIError on API errors.

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand, gene_id,
        biotype, assembly. Empty on API error unless raise_on_error=True.

    Raises:
        UCSCAPIError: If raise_on_error=True and the API fails.
    """
    empty = pd.DataFrame(columns=_gene_columns())
    try:
        rows = _fetch_track(ucsc_genome, chrom, start, end)
    except UCSCAPIError:
        if raise_on_error:
            raise
        return empty

    assembly = _assembly_for(ucsc_genome)
    chrom_str = normalize_chrom(chrom)

    genes: dict[str, dict] = {}
    for row in rows:
        symbol = row.get("name2") or row.get("name")
        if not symbol:
            continue
        coding = str(row.get("name", "")).startswith(_CODING_PREFIXES)
        gene = genes.get(symbol)
        if gene is None:
            genes[symbol] = {
                "chr": chrom_str,
                "start": int(row["txStart"]) + 1,
                "end": int(row["txEnd"]),
                "gene_name": symbol,
                "strand": row.get("strand", "+"),
                "gene_id": symbol,
                "biotype": "protein_coding" if coding else "non_coding",
                "assembly": assembly,
            }
            continue
        gene["start"] = min(gene["start"], int(row["txStart"]) + 1)
        gene["end"] = max(gene["end"], int(row["txEnd"]))
        if coding:
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
    empty = pd.DataFrame(columns=_exon_columns())
    try:
        rows = _fetch_track(ucsc_genome, chrom, start, end)
    except UCSCAPIError:
        if raise_on_error:
            raise
        return empty

    assembly = _assembly_for(ucsc_genome)
    chrom_str = normalize_chrom(chrom)

    records = []
    for row in rows:
        transcript = str(row.get("name", ""))
        starts = str(row.get("exonStarts", "")).strip(",")
        ends = str(row.get("exonEnds", "")).strip(",")
        if not starts or not ends:
            continue
        for index, (exon_start, exon_end) in enumerate(
            zip(starts.split(","), ends.split(","))
        ):
            records.append(
                {
                    "chr": chrom_str,
                    "start": int(exon_start) + 1,
                    "end": int(exon_end),
                    "gene_name": row.get("name2", ""),
                    "exon_id": f"{transcript}_exon{index + 1}",
                    "transcript_id": transcript,
                    "assembly": assembly,
                }
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
    cache_dir=None,
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
        genes_df = pd.DataFrame(columns=_gene_columns())

    if use_cache and not fetch_failed:
        save_genes(genes_df, cache_dir, ucsc_genome, chrom_str, start, end)

    return (genes_df, exons()) if include_exons else genes_df


def clear_ucsc_cache(cache_dir=None, ucsc_genome: str | None = None) -> int:
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

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

import pandas as pd

from ._gene_source import (
    EXON_COLUMNS,
    GENE_COLUMNS,
    GeneAnnotations,
    GeneSource,
    empty_frame,
)
from ._http import request_json
from .exceptions import UCSCAPIError, ValidationError
from .genome_build import GenomeBuild, resolve_build, ucsc_chrom
from .logging import logger
from .utils import normalize_chrom

UCSC_REST_URL = "https://api.genome.ucsc.edu"
UCSC_GENE_TRACK = "ncbiRefSeq"

# RefSeq accession prefixes for transcripts that code for protein.
_CODING_PREFIXES = ("NM_", "XM_")


def _ucsc_build(build: str | GenomeBuild) -> GenomeBuild:
    """Resolve a build that UCSC serves gene annotations for."""
    record = resolve_build(build)
    if record is None or record.ucsc_genome is None:
        raise ValidationError(f"UCSC serves no gene annotations for build {build!r}")
    return record


def _fetch_track(
    build: GenomeBuild, chrom: str | int, start: int, end: int
) -> list[dict]:
    """Fetch raw ncbiRefSeq transcript rows overlapping a region."""
    ucsc_name = ucsc_chrom(chrom, build)
    payload = request_json(
        f"{UCSC_REST_URL}/getData/track",
        {
            "genome": build.ucsc_genome,
            "track": UCSC_GENE_TRACK,
            "chrom": ucsc_name,
            "start": start,
            "end": end,
        },
        error_cls=UCSCAPIError,
        service="UCSC",
    )
    rows = payload.get(UCSC_GENE_TRACK, [])
    # UCSC returns a chrom-keyed dict when the request spans the whole genome.
    if isinstance(rows, dict):
        rows = rows.get(ucsc_name, [])
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


def _genes_from_rows(
    rows: list[dict], build: GenomeBuild, chrom_str: str, biotype: str
) -> pd.DataFrame:
    """Collapse transcript rows into one row per gene symbol.

    ncbiRefSeq is a transcript-level track, so the many transcripts sharing a
    symbol become one row spanning the widest of them.
    """
    assembly = build.assembly_name

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
        logger.debug(f"No genes found in {build.ucsc_genome} {chrom_str}")
        return empty_frame(GENE_COLUMNS)

    logger.debug(f"Fetched {len(records)} genes from UCSC {build.ucsc_genome}")
    return pd.DataFrame(records)


def _exons_from_rows(
    rows: list[dict], build: GenomeBuild, chrom_str: str
) -> pd.DataFrame:
    """Expand every transcript row into one row per exon."""
    assembly = build.assembly_name

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
        logger.debug(f"No exons found in {build.ucsc_genome} {chrom_str}")
        return empty_frame(EXON_COLUMNS)

    logger.debug(f"Fetched {len(records)} exons from UCSC {build.ucsc_genome}")
    return pd.DataFrame(records)


def fetch_track_frames(
    build: str | GenomeBuild,
    chrom: str | int,
    start: int,
    end: int,
    biotype: str = "protein_coding",
) -> GeneAnnotations:
    """Fetch the genes and exons for a region from one track request.

    ncbiRefSeq rows already carry ``exonStarts``/``exonEnds``, so the exons
    cost nothing beyond the genes.

    Args:
        build: Build UCSC serves genes for, as a record or a name such as
            ``"canFam3"`` or ``"CanFam3.1"``.
        chrom: Chromosome name or number. Codes UCSC spells differently,
            such as canine 39, are asked for by their UCSC name; the rows keep
            the caller's name.
        start: Region start position (1-based).
        end: Region end position (1-based).
        biotype: Gene biotype filter. ncbiRefSeq carries no finer biotypes
            than the accession prefix, so the only meaningful values are
            ``"protein_coding"`` (at least one NM_/XM_ transcript),
            ``"non_coding"`` (none), and None or "" to keep everything; any
            other value matches nothing.

    Returns:
        The region's annotations. UCSC coordinates are 0-based half-open
        and are converted to the 1-based inclusive convention Ensembl and the
        rest of pyLocusZoom use.

    Raises:
        ValidationError: If UCSC serves no gene annotations for the build.
        UCSCAPIError: If the API fails.
    """
    record = _ucsc_build(build)
    chrom_str = normalize_chrom(chrom)
    rows = _fetch_track(record, chrom_str, start, end)

    return GeneAnnotations(
        _genes_from_rows(rows, record, chrom_str, biotype),
        _exons_from_rows(rows, record, chrom_str),
    )


def ucsc_source(build: str | GenomeBuild) -> GeneSource:
    """Build the GeneSource for one build UCSC serves genes for."""
    record = _ucsc_build(build)
    return GeneSource(
        name="ucsc",
        cache_species=record.ucsc_genome,
        build_token="",
        fetch=lambda chrom, start, end: fetch_track_frames(record, chrom, start, end),
    )

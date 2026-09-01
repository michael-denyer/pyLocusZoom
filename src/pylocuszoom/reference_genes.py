# src/pylocuszoom/reference_genes.py
"""Pick the gene-annotation source that can serve the caller's genome build.

Ensembl serves exactly one reference assembly per species and answers a request
naming any other with that same assembly and an HTTP 200, so it cannot supply
genes for a build it has retired and will not say so. Release 116 was the last
on the legacy REST platform and the archive REST hosts redirect to a help page,
so CanFam3.1, CanFam4 and FelCat9 have no Ensembl source at any URL. UCSC hosts
all three.

``UCSC_BUILDS`` is the whole policy: a build listed there is fetched from UCSC,
anything else from Ensembl. Both sources return the same columns, including
``assembly``, so callers do not branch on which one answered.

A ``GeneSource`` is everything the fetch-and-cache orchestration needs to know
about a source: where it caches, which error it raises, and how to ask it for
genes and exons. ``get_genes_for_build`` holds the one copy of that
orchestration, so the cache policy and the error translation exist once rather
than once per source.
"""

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import pandas as pd

from ._gene_cache import cache_root, clear_cache, load_genes, save_genes
from .exceptions import EnsemblAPIError, ReferenceAPIError, UCSCAPIError
from .logging import logger
from .utils import assembly_token, normalize_chrom

# Build token (from utils.assembly_token) -> UCSC genome name.
# Only builds Ensembl cannot serve belong here; everything else stays on
# Ensembl, which covers far more species than UCSC does.
UCSC_BUILDS: dict[str, str] = {
    "canfam3": "canFam3",
    "canfam4": "canFam4",
    "felcat9": "felCat9",
}

# The frame schema both sources are contractually required to produce. The
# cache CSV round-trip depends on it: a column-less empty frame serialises to a
# one-byte file that pd.read_csv cannot parse back.
GENE_COLUMNS = (
    "chr",
    "start",
    "end",
    "gene_name",
    "strand",
    "gene_id",
    "biotype",
    "assembly",
)
EXON_COLUMNS = (
    "chr",
    "start",
    "end",
    "gene_name",
    "exon_id",
    "transcript_id",
    "assembly",
)

# What a fetch is being asked for. "both" exists so UCSC can serve genes and
# exons from one track request; Ensembl needs a request per feature type either
# way, and the cache-hit path asks for "exons" alone.
FetchWhat = Literal["genes", "exons", "both"]

# Fetch (genes, exons) for (chrom, start, end, what). The frame that was not
# asked for comes back empty.
FetchFn = Callable[[str, int, int, FetchWhat], "tuple[pd.DataFrame, pd.DataFrame]"]


@dataclass(frozen=True)
class GeneSource:
    """One gene-annotation source, parametrising the shared orchestration.

    Attributes:
        name: Cache leaf directory, e.g. ``"ensembl"`` or ``"ucsc"``.
        error_cls: Error the fetchers raise when the service fails.
        cache_species: Key the cache is partitioned by, e.g. the Ensembl
            species name or the UCSC genome name.
        build_token: Extra cache-key component. Ensembl serves every build of a
            species from one URL, so its entries have to be keyed by build;
            a UCSC genome already names one build, so it needs none.
        fetch: Fetch (genes, exons) for (chrom, start, end, what). Always
            raises ``error_cls`` on failure rather than returning an empty
            frame. The frame that was not asked for comes back empty.
        on_cache_hit: Inspect a frame reloaded from cache. Ensembl uses it to
            repeat its assembly-mismatch warning in a later session.
    """

    name: str
    error_cls: type[ReferenceAPIError]
    cache_species: str
    build_token: str
    fetch: FetchFn
    on_cache_hit: Callable[[pd.DataFrame], None] = lambda cached: None


def ucsc_genome_for_build(genome_build: str | None) -> str | None:
    """Return the UCSC genome serving this build, or None to use Ensembl."""
    if not genome_build:
        return None
    return UCSC_BUILDS.get(assembly_token(genome_build))


def ucsc_source(ucsc_genome: str) -> GeneSource:
    """Build the GeneSource for one UCSC genome."""
    from .ucsc import fetch_track_frames

    return GeneSource(
        name="ucsc",
        error_cls=UCSCAPIError,
        cache_species=ucsc_genome,
        build_token="",
        fetch=lambda chrom, start, end, what: fetch_track_frames(
            ucsc_genome, chrom, start, end, what
        ),
    )


def ensembl_source(species: str, genome_build: str | None = None) -> GeneSource:
    """Build the GeneSource for one species on Ensembl's current assembly."""
    from .ensembl import (
        fetch_overlap_frames,
        get_ensembl_species_name,
        warn_on_cached_assembly,
    )

    ensembl_species = get_ensembl_species_name(species)
    return GeneSource(
        name="ensembl",
        error_cls=EnsemblAPIError,
        cache_species=ensembl_species,
        build_token=assembly_token(genome_build or ""),
        fetch=lambda chrom, start, end, what: fetch_overlap_frames(
            species, chrom, start, end, what, genome_build=genome_build
        ),
        on_cache_hit=lambda cached: warn_on_cached_assembly(
            cached, genome_build, ensembl_species
        ),
    )


def source_for(species: str, genome_build: str | None = None) -> GeneSource:
    """Pick the GeneSource that can serve this species and genome build.

    Args:
        species: Species name or alias, used when Ensembl serves the build.
        genome_build: Build the caller's data is in. A build only UCSC serves
            selects the UCSC source; anything else stays on Ensembl.

    Returns:
        The source to fetch from, already bound to the species or genome.
    """
    ucsc_genome = ucsc_genome_for_build(genome_build)
    if ucsc_genome is None:
        return ensembl_source(species, genome_build)

    logger.debug(
        f"Genome build {genome_build!r} is not served by Ensembl; "
        f"fetching genes from UCSC {ucsc_genome}"
    )
    return ucsc_source(ucsc_genome)


def get_genes_for_build(
    species: str,
    chrom: str | int,
    start: int,
    end: int,
    genome_build: str | None = None,
    cache_dir: Path | None = None,
    use_cache: bool = True,
    include_exons: bool = False,
    raise_on_error: bool = False,
    source: GeneSource | None = None,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Get gene annotations for a region in the caller's genome build.

    Checks the disk cache first and fetches from the serving API on a miss.
    An empty region is cached so gene-sparse regions stop re-requesting; a
    failed fetch is never cached, because on reload it would be
    indistinguishable from an empty region and would permanently hide the
    region's genes.

    Args:
        species: Species name or alias, used when Ensembl serves the build.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        genome_build: Build the caller's data is in. Selects the source; None
            means Ensembl's current assembly for the species.
        cache_dir: Cache directory (the source's default if None).
        use_cache: Whether to use the disk cache.
        include_exons: If True, return a (genes_df, exons_df) tuple.
        raise_on_error: If True, raise on API errors instead of returning empty.
        source: Source to fetch from. Derived from species and genome_build if
            None; the public per-source entry points pass their own.

    Returns:
        Gene annotations in the requested build, or a (genes_df, exons_df)
        tuple when include_exons. Both carry an ``assembly`` column.

    Raises:
        ValidationError: If an Ensembl-served build asks for a region over
            5Mb. UCSC imposes no region limit, so UCSC-served builds don't.
        ReferenceAPIError: If raise_on_error=True and the serving API fails
            (EnsemblAPIError or UCSCAPIError, depending on the source).
    """
    if source is None:
        source = source_for(species, genome_build)
    if cache_dir is None:
        cache_dir = cache_root(source.name)

    chrom_str = normalize_chrom(chrom)
    empty_genes = pd.DataFrame(columns=list(GENE_COLUMNS))
    empty_exons = pd.DataFrame(columns=list(EXON_COLUMNS))

    def fetch(what: FetchWhat) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
        try:
            genes_df, exons_df = source.fetch(chrom_str, start, end, what)
        except source.error_cls:
            if raise_on_error:
                raise
            return empty_genes, empty_exons, True
        return genes_df, exons_df, False

    if use_cache:
        cached = load_genes(
            cache_dir,
            source.cache_species,
            chrom_str,
            start,
            end,
            build_token=source.build_token,
        )
        if cached is not None:
            source.on_cache_hit(cached)
            if not include_exons:
                return cached
            return cached, fetch("exons")[1]

    genes_df, exons_df, fetch_failed = fetch("both" if include_exons else "genes")

    if use_cache and not fetch_failed:
        save_genes(
            genes_df,
            cache_dir,
            source.cache_species,
            chrom_str,
            start,
            end,
            build_token=source.build_token,
        )

    return (genes_df, exons_df) if include_exons else genes_df


def clear_gene_cache(
    source_name: str,
    cache_dir: Path | None = None,
    cache_species: str | None = None,
) -> int:
    """Clear one source's cached gene files.

    Args:
        source_name: Cache leaf directory, e.g. ``"ensembl"`` or ``"ucsc"``.
        cache_dir: Cache directory (the source's default if None).
        cache_species: If given, only this species or genome subdirectory is
            cleared.

    Returns:
        Number of files deleted.
    """
    if cache_dir is None:
        cache_dir = cache_root(source_name)
    deleted = clear_cache(cache_dir, cache_species)
    logger.info(f"Cleared {deleted} cached {source_name} files from {cache_dir}")
    return deleted

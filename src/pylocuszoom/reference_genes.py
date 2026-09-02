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
about a source: where it caches and how to ask it for genes and exons.
``get_genes_for_build`` holds the one copy of that orchestration, so the cache
policy exists once rather than once per source.
"""

from pathlib import Path

from ._gene_cache import cache_root, clear_cache, load_annotations, save_annotations
from ._gene_source import GeneAnnotations, GeneSource
from .ensembl import ensembl_source
from .logging import logger
from .ucsc import ucsc_source
from .utils import assembly_token, normalize_chrom

# Build token (from utils.assembly_token) -> UCSC genome name.
# Only builds Ensembl cannot serve belong here; everything else stays on
# Ensembl, which covers far more species than UCSC does.
UCSC_BUILDS: dict[str, str] = {
    "canfam3": "canFam3",
    "canfam4": "canFam4",
    "felcat9": "felCat9",
}


def ucsc_genome_for_build(genome_build: str | None) -> str | None:
    """Return the UCSC genome serving this build, or None to use Ensembl."""
    if not genome_build:
        return None
    return UCSC_BUILDS.get(assembly_token(genome_build))


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
    source: GeneSource,
    chrom: str | int,
    start: int,
    end: int,
    *,
    cache_dir: Path | None = None,
    use_cache: bool = True,
) -> GeneAnnotations:
    """Get the gene annotations for a region from one source.

    Checks the disk cache first and fetches from the serving API on a miss.
    An empty region is cached so gene-sparse regions stop re-requesting; a
    failed fetch is never cached, because on reload it would be
    indistinguishable from an empty region and would permanently hide the
    region's genes.

    Args:
        source: Source to fetch from, from ``source_for``.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        cache_dir: Cache directory (the source's default if None).
        use_cache: Whether to use the disk cache.

    Returns:
        The region's genes and exons, both carrying an ``assembly`` column.

    Raises:
        ValidationError: If an Ensembl-served build asks for a region over
            5Mb. UCSC imposes no region limit, so UCSC-served builds don't.
        ReferenceAPIError: If the serving API fails (EnsemblAPIError or
            UCSCAPIError, depending on the source).
    """
    if cache_dir is None:
        cache_dir = cache_root(source.name)

    chrom_str = normalize_chrom(chrom)
    key = (cache_dir, source.cache_species, chrom_str, start, end)

    if use_cache:
        cached = load_annotations(*key, build_token=source.build_token)
        if cached is not None:
            source.on_cache_hit(cached.genes)
            return cached

    annotations = source.fetch(chrom_str, start, end)

    if use_cache:
        save_annotations(annotations, *key, build_token=source.build_token)

    return annotations


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

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
"""

from pathlib import Path

import pandas as pd

from .ensembl import get_genes_for_region
from .logging import logger
from .ucsc import get_genes_for_region_ucsc
from .utils import assembly_token

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
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Get gene annotations for a region in the caller's genome build.

    Args:
        species: Species name or alias, used when Ensembl serves the build.
        chrom: Chromosome name or number.
        start: Region start position (1-based).
        end: Region end position (1-based).
        genome_build: Build the caller's data is in. Selects the source; None
            means Ensembl's current assembly for the species.
        cache_dir: Cache directory (each source uses its own default if None).
        use_cache: Whether to use the disk cache.
        include_exons: If True, return a (genes_df, exons_df) tuple.
        raise_on_error: If True, raise on API errors instead of returning empty.

    Returns:
        Gene annotations in the requested build, or a (genes_df, exons_df)
        tuple when include_exons. Both carry an ``assembly`` column.

    Raises:
        ValidationError: If an Ensembl-served build asks for a region over
            5Mb. UCSC imposes no region limit, so UCSC-served builds don't.
        ReferenceAPIError: If raise_on_error=True and the serving API fails
            (EnsemblAPIError or UCSCAPIError, depending on the source).
    """
    ucsc_genome = ucsc_genome_for_build(genome_build)
    if ucsc_genome is not None:
        logger.debug(
            f"Genome build {genome_build!r} is not served by Ensembl; "
            f"fetching genes from UCSC {ucsc_genome}"
        )
        return get_genes_for_region_ucsc(
            ucsc_genome,
            chrom,
            start,
            end,
            cache_dir=cache_dir,
            use_cache=use_cache,
            include_exons=include_exons,
            raise_on_error=raise_on_error,
        )

    return get_genes_for_region(
        species,
        chrom,
        start,
        end,
        cache_dir=cache_dir,
        use_cache=use_cache,
        include_exons=include_exons,
        raise_on_error=raise_on_error,
        genome_build=genome_build,
    )

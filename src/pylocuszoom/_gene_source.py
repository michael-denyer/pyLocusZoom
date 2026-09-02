# src/pylocuszoom/_gene_source.py
"""The vocabulary both gene sources speak, in a leaf neither of them imports back.

``ensembl.py`` and ``ucsc.py`` each export one ``GeneSource`` constructor, and
``reference_genes.py`` picks between them. Holding the frame schema and the
value type here rather than in ``reference_genes.py`` is what lets the routing
module import its two sources at the top of the file instead of inside its
functions.
"""

from collections.abc import Callable
from dataclasses import dataclass

import pandas as pd

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

# Fetch (genes, exons) for (chrom, start, end). Both sources answer from one
# request, so neither frame costs anything the other did not.
FetchFn = Callable[[str, int, int], "tuple[pd.DataFrame, pd.DataFrame]"]


def empty_frame(columns: tuple[str, ...]) -> pd.DataFrame:
    """Build the empty frame for a schema, columns and all."""
    return pd.DataFrame(columns=list(columns))


@dataclass(frozen=True)
class GeneSource:
    """One gene-annotation source, parametrising the shared orchestration.

    Attributes:
        name: Cache leaf directory, e.g. ``"ensembl"`` or ``"ucsc"``.
        cache_species: Key the cache is partitioned by, e.g. the Ensembl
            species name or the UCSC genome name.
        build_token: Extra cache-key component. Ensembl serves every build of a
            species from one URL, so its entries have to be keyed by build;
            a UCSC genome already names one build, so it needs none.
        fetch: Fetch (genes, exons) for (chrom, start, end). Always raises
            ``ReferenceAPIError`` on failure rather than returning an empty
            frame.
        on_cache_hit: Inspect a frame reloaded from cache. Ensembl uses it to
            repeat its assembly-mismatch warning in a later session.
    """

    name: str
    cache_species: str
    build_token: str
    fetch: FetchFn
    on_cache_hit: Callable[[pd.DataFrame], None] = lambda cached: None

"""File format loaders for common GWAS, eQTL, and fine-mapping outputs.

Convenience functions to load data from standard file formats into
DataFrames ready for use with LocusZoomPlotter. One module per family, each
declaring its formats as :class:`~._engine.LoaderSpec` constants that the
shared ``_load_tabular`` engine runs.

GWAS formats (:mod:`~.gwas`):
- PLINK (.assoc, .assoc.linear, .assoc.logistic, .qassoc)
- REGENIE (.regenie)
- BOLT-LMM (.stats)
- GEMMA (.assoc.txt)
- SAIGE (.txt)
- Generic TSV/CSV

eQTL formats (:mod:`~.eqtl`):
- GTEx significant pairs format
- eQTL Catalogue format
- MatrixEQTL output

Fine-mapping formats (:mod:`~.finemapping`):
- SuSiE (susieR output)
- FINEMAP (.snp output)
- CAVIAR (.set output)

Gene annotation formats (:mod:`~.annotation`):
- GTF/GFF3
- BED (4-column: chr, start, end, name)
"""

from .annotation import load_bed, load_ensembl_genes, load_gtf
from .eqtl import load_eqtl_catalogue, load_gtex_eqtl, load_matrixeqtl
from .finemapping import load_caviar, load_finemap, load_polyfun, load_susie
from .gwas import (
    load_bolt_lmm,
    load_gemma,
    load_gwas,
    load_gwas_catalog,
    load_plink_assoc,
    load_regenie,
    load_saige,
)

__all__ = [
    "load_bed",
    "load_bolt_lmm",
    "load_caviar",
    "load_ensembl_genes",
    "load_eqtl_catalogue",
    "load_finemap",
    "load_gemma",
    "load_gtex_eqtl",
    "load_gtf",
    "load_gwas",
    "load_gwas_catalog",
    "load_matrixeqtl",
    "load_plink_assoc",
    "load_polyfun",
    "load_regenie",
    "load_saige",
    "load_susie",
]

"""Gene annotation loaders: GTF/GFF3, BED, and Ensembl BioMart exports.

Each returns a DataFrame with columns chr, start, end, gene_name, and strand
where the format carries one.

Raises:
    LoaderValidationError: If the format's columns cannot be mapped or the
        mapped values fail the format's schema.
"""

from dataclasses import replace
from pathlib import Path
from typing import Union

import pandas as pd

from ..logging import logger
from ..schemas import Family, Tier, spec
from ..utils import normalize_chrom_series
from ..validation import check
from ._engine import LoaderSpec, _load_tabular


def load_gtf(
    filepath: Union[str, Path],
    feature_type: str = "gene",
) -> pd.DataFrame:
    """Load gene annotations from GTF/GFF3 file.

    Args:
        filepath: Path to GTF or GFF3 file (can be gzipped).
        feature_type: Feature type to extract ("gene", "exon", "transcript").
            Default "gene".

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand.

    Example:
        >>> genes_df = load_gtf("genes.gtf", feature_type="gene")
        >>> exons_df = load_gtf("genes.gtf", feature_type="exon")
    """
    # GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
    df = pd.read_csv(
        filepath,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
    )

    # Filter to requested feature type
    df = df[df["feature"] == feature_type].copy()

    # Parse gene_name from attributes
    def extract_gene_name(attrs: str) -> str:
        """Extract gene_name or gene_id from GTF attributes."""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                # gene_name "BRCA1" or gene_name=BRCA1
                return attr.split('"')[1] if '"' in attr else attr.split("=")[1]
            if attr.startswith("gene_id"):
                return attr.split('"')[1] if '"' in attr else attr.split("=")[1]
        return ""

    df["gene_name"] = df["attributes"].apply(extract_gene_name)

    # Clean chromosome names
    df["chr"] = normalize_chrom_series(df["chr"])

    # Select and return relevant columns
    result = df[["chr", "start", "end", "gene_name", "strand"]].copy()
    logger.debug(f"Loaded {len(result)} {feature_type} features from GTF")
    check(result, spec(Family.GENES, Tier.LOAD))
    return result


# Standard BED column names, in order, up to BED12.
_BED_COLUMNS = (
    "chr",
    "start",
    "end",
    "gene_name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
)


def _bed_names(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Name a headerless BED frame, generically past BED12."""
    if not all(isinstance(c, int) for c in df.columns):
        return df
    n_cols = len(df.columns)
    names = list(_BED_COLUMNS[:n_cols])
    names += [f"col{i}" for i in range(len(_BED_COLUMNS), n_cols)]
    df.columns = names
    return df


_BED_SPEC = LoaderSpec(
    log_fmt="Loaded {n} features from BED",
    read={"sep": "\t"},
    col_map={
        "chrom": "chr",
        "chromStart": "start",
        "chromEnd": "end",
        "name": "gene_name",
    },
    transform=_bed_names,
    clean_chrom=True,
    schema=lambda out_cols: spec(Family.GENES, Tier.LOAD),
)


def load_bed(
    filepath: Union[str, Path],
    has_header: bool = False,
) -> pd.DataFrame:
    """Load gene annotations from BED file.

    Supports BED4+ format (chr, start, end, name, ...).

    Args:
        filepath: Path to BED file.
        has_header: Whether file has header row. Default False.

    Returns:
        DataFrame with columns: chr, start, end, gene_name.

    Example:
        >>> genes_df = load_bed("genes.bed")
    """
    spec = replace(
        _BED_SPEC, read={**_BED_SPEC.read, "header": 0 if has_header else None}
    )
    return _load_tabular(filepath, spec)


def _ensembl_strand(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Convert Ensembl 1/-1 strand encoding to +/- symbols."""
    if "strand" in df.columns:
        df["strand"] = df["strand"].map({1: "+", -1: "-", "+": "+", "-": "-"})
    return df


_ENSEMBL_SPEC = LoaderSpec(
    log_fmt="Loaded {n} genes from Ensembl export",
    read={"sep": "\t"},
    col_map={
        "Chromosome/scaffold name": "chr",
        "Gene start (bp)": "start",
        "Gene end (bp)": "end",
        "Gene name": "gene_name",
        "Strand": "strand",
        # Alternative column names
        "chromosome_name": "chr",
        "start_position": "start",
        "end_position": "end",
        "external_gene_name": "gene_name",
    },
    transform=_ensembl_strand,
    schema=lambda out_cols: spec(Family.GENES, Tier.LOAD),
)


def load_ensembl_genes(
    filepath: Union[str, Path],
) -> pd.DataFrame:
    """Load Ensembl BioMart gene export.

    Args:
        filepath: Path to BioMart export file (TSV).

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand.
    """
    return _load_tabular(filepath, _ENSEMBL_SPEC)

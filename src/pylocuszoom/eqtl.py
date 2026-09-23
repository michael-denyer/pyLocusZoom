"""eQTL data handling and validation for pyLocusZoom.

Provides utilities for loading, validating, and preparing expression
quantitative trait loci (eQTL) data for overlay on regional plots.
"""

from typing import List, Optional

import numpy as np
import pandas as pd

from ._data import prepare_pvalue_data
from .exceptions import EQTLValidationError
from .logging import logger
from .schemas import Canonical, eqtl_plot_spec
from .utils import filter_by_region, normalize_chrom, normalize_chrom_series
from .validation import ColumnSpec, check


def filter_eqtl_by_gene(
    df: pd.DataFrame,
    gene: str,
    gene_col: str = "gene",
) -> pd.DataFrame:
    """Filter eQTL data to a specific target gene.

    Args:
        df: eQTL DataFrame.
        gene: Target gene name to filter for.
        gene_col: Column containing gene names.

    Returns:
        Filtered DataFrame containing only eQTLs for the target gene.

    Raises:
        EQTLValidationError: If gene column doesn't exist.
    """
    if gene_col not in df.columns:
        raise EQTLValidationError(
            f"Cannot filter by gene: column '{gene_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    filtered = df[df[gene_col] == gene].copy()
    logger.debug(f"Filtered eQTL data to {len(filtered)} variants for gene {gene}")
    return filtered


def filter_eqtl_by_region(
    df: pd.DataFrame,
    chrom: int,
    start: int,
    end: int,
    pos_col: str = Canonical.POS,
    chrom_col: Optional[str] = Canonical.CHROM,
) -> pd.DataFrame:
    """Filter eQTL data to a genomic region.

    Args:
        df: eQTL DataFrame.
        chrom: Chromosome number.
        start: Start position.
        end: End position.
        pos_col: Column name for position.
        chrom_col: Column name for chromosome, or None to filter by position
            only.

    Returns:
        Filtered DataFrame containing only eQTLs in the region.
    """
    filtered = filter_by_region(
        df,
        region=(chrom, start, end),
        chrom_col=chrom_col,
        pos_col=pos_col,
    )
    logger.debug(
        f"Filtered eQTL data to {len(filtered)} variants in region chr{chrom}:{start}-{end}"
    )
    return filtered


def prepare_eqtl_for_plotting(
    df: pd.DataFrame,
    pos_col: str = Canonical.POS,
    p_col: str = Canonical.P,
    gene: Optional[str] = None,
    chrom: Optional[int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    chrom_col: Optional[str] = Canonical.CHROM,
) -> pd.DataFrame:
    """Prepare eQTL data for plotting.

    Validates, filters, and adds computed columns needed for plotting.

    Args:
        df: Raw eQTL DataFrame.
        pos_col: Column name for position.
        p_col: Column name for p-value.
        gene: Optional gene to filter for.
        chrom: Optional chromosome for region filtering.
        start: Optional start position for region filtering.
        end: Optional end position for region filtering.
        chrom_col: Chromosome column for region filtering, or None to filter
            by position only.

    Returns:
        Prepared DataFrame with neglog10p column added.
    """
    check(df, eqtl_plot_spec(pos_col, p_col))

    result = df.copy()

    # Filter by gene if specified
    if gene:
        result = filter_eqtl_by_gene(result, gene)

    # Filter by region if specified
    if chrom is not None and start is not None and end is not None:
        result = filter_eqtl_by_region(
            result, chrom, start, end, pos_col=pos_col, chrom_col=chrom_col
        )

    return prepare_pvalue_data(result, p_col, "eqtl")


def get_eqtl_genes(df: pd.DataFrame, gene_col: str = "gene") -> List[str]:
    """Get list of unique genes in eQTL data.

    Args:
        df: eQTL DataFrame.
        gene_col: Column containing gene names.

    Returns:
        Sorted list of unique gene names.
    """
    if gene_col not in df.columns:
        return []
    return sorted(df[gene_col].dropna().unique().tolist())


def _overlap_coordinates(
    df: pd.DataFrame,
    pos_col: str,
    p_col: str,
    chrom_col: str,
    common_chrom: int | str | None,
) -> pd.DataFrame:
    """Resolve each input to coordinate and p-value roles before joining."""
    check(
        df,
        ColumnSpec(
            name="eQTL DataFrame",
            required=(pos_col, p_col),
            numeric=(p_col,),
            error_class=EQTLValidationError,
        ),
    )
    positions = pd.to_numeric(df[pos_col], errors="coerce")
    if (
        df[pos_col].map(lambda value: isinstance(value, (bool, np.bool_))).any()
        or pd.api.types.is_complex_dtype(positions)
        or positions.isna().any()
        or not np.isfinite(positions).all()
        or not positions.between(1, 2**63, inclusive="left").all()
        or positions.mod(1).ne(0).any()
    ):
        raise EQTLValidationError(
            f"Absolute positions in {pos_col!r} must be finite positive integers "
            "below 2**63"
        )
    positions = positions.astype("int64")
    if chrom_col in df.columns:
        if df[chrom_col].isna().any():
            raise EQTLValidationError("Chromosome coordinates must not be null")
        chromosomes = normalize_chrom_series(df[chrom_col])
        if common_chrom is not None and not chromosomes.eq(common_chrom).all():
            raise EQTLValidationError(
                "Input chromosome coordinates disagree with common_chrom"
            )
    elif common_chrom is not None:
        chromosomes = common_chrom
    else:
        raise EQTLValidationError(
            f"Missing chromosome column {chrom_col!r}. For inputs already scoped "
            "to the same chromosome, supply common_chrom explicitly."
        )
    return pd.DataFrame({"chr": chromosomes, "pos": positions, "p_value": df[p_col]})


def calculate_colocalization_overlap(
    gwas_df: pd.DataFrame,
    eqtl_df: pd.DataFrame,
    gwas_pos_col: str = Canonical.POS,
    eqtl_pos_col: str = Canonical.POS,
    gwas_p_col: str = Canonical.P,
    eqtl_p_col: str = Canonical.P,
    p_threshold: float = 1e-5,
    *,
    gwas_chrom_col: str = Canonical.CHROM,
    eqtl_chrom_col: str = Canonical.CHROM,
    common_chrom: int | str | None = None,
) -> pd.DataFrame:
    """Find significant GWAS/eQTL overlaps by chromosome and absolute position.

    This is coordinate overlap, without allele matching or harmonization.
    For formal colocalization use dedicated tools such as coloc or eCAVIAR.

    Args:
        gwas_df: GWAS results DataFrame.
        eqtl_df: eQTL results DataFrame.
        gwas_pos_col: Absolute position column in GWAS data.
        eqtl_pos_col: Absolute position column in eQTL data.
        gwas_p_col: P-value column in GWAS data.
        eqtl_p_col: P-value column in eQTL data.
        p_threshold: P-value threshold for significance.
        gwas_chrom_col: Chromosome column in GWAS data.
        eqtl_chrom_col: Chromosome column in eQTL data.
        common_chrom: Explicit shared chromosome for position-only inputs.
            Any chromosome columns present must agree with this value.

    Returns:
        DataFrame with chr, pos, p_value_gwas and p_value_eqtl columns.
        Multiple variants at one coordinate are retained as multiple matches;
        this does not establish allele-level identity.

    Raises:
        EQTLValidationError: If positions are not finite positive integers,
            or chromosomes are missing or contradict the declared common chromosome.
    """
    chromosome = normalize_chrom(common_chrom) if common_chrom is not None else None
    gwas = _overlap_coordinates(
        gwas_df, gwas_pos_col, gwas_p_col, gwas_chrom_col, chromosome
    )
    eqtl = _overlap_coordinates(
        eqtl_df, eqtl_pos_col, eqtl_p_col, eqtl_chrom_col, chromosome
    )
    overlap = gwas[gwas["p_value"] < p_threshold].merge(
        eqtl[eqtl["p_value"] < p_threshold],
        on=["chr", "pos"],
        how="inner",
        suffixes=("_gwas", "_eqtl"),
    )
    logger.info(
        f"Found {len(overlap)} coordinate matches significant in both GWAS and eQTL "
        f"(p < {p_threshold})"
    )
    return overlap

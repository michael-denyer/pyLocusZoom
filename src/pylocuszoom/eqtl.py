"""eQTL data handling and validation for pyLocusZoom.

Provides utilities for loading, validating, and preparing expression
quantitative trait loci (eQTL) data for overlay on regional plots.
"""

from typing import List, Optional

import numpy as np
import pandas as pd

from .logging import logger

REQUIRED_EQTL_COLS = ["pos", "p_value"]
OPTIONAL_EQTL_COLS = ["gene", "effect_size", "rs", "se"]


class EQTLValidationError(ValueError):
    """Raised when eQTL DataFrame validation fails."""

    pass


def validate_eqtl_df(
    df: pd.DataFrame,
    pos_col: str = "pos",
    p_col: str = "p_value",
) -> None:
    """Validate eQTL DataFrame has required columns.

    Args:
        df: eQTL DataFrame to validate.
        pos_col: Column name for genomic position.
        p_col: Column name for p-value.

    Raises:
        EQTLValidationError: If required columns are missing.
    """
    missing = []
    if pos_col not in df.columns:
        missing.append(pos_col)
    if p_col not in df.columns:
        missing.append(p_col)

    if missing:
        raise EQTLValidationError(
            f"eQTL DataFrame missing required columns: {missing}. "
            f"Required: {pos_col} (position), {p_col} (p-value)"
        )


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
    pos_col: str = "pos",
    chrom_col: Optional[str] = "chr",
) -> pd.DataFrame:
    """Filter eQTL data to a genomic region.

    Args:
        df: eQTL DataFrame.
        chrom: Chromosome number.
        start: Start position.
        end: End position.
        pos_col: Column name for position.
        chrom_col: Column name for chromosome (if present).

    Returns:
        Filtered DataFrame containing only eQTLs in the region.
    """
    mask = (df[pos_col] >= start) & (df[pos_col] <= end)

    # Filter by chromosome if column exists
    if chrom_col and chrom_col in df.columns:
        chrom_str = str(chrom).replace("chr", "")
        df_chrom = df[chrom_col].astype(str).str.replace("chr", "", regex=False)
        mask = mask & (df_chrom == chrom_str)

    filtered = df[mask].copy()
    logger.debug(f"Filtered eQTL data to {len(filtered)} variants in region chr{chrom}:{start}-{end}")
    return filtered


def prepare_eqtl_for_plotting(
    df: pd.DataFrame,
    pos_col: str = "pos",
    p_col: str = "p_value",
    gene: Optional[str] = None,
    chrom: Optional[int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
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

    Returns:
        Prepared DataFrame with neglog10p column added.
    """
    validate_eqtl_df(df, pos_col=pos_col, p_col=p_col)

    result = df.copy()

    # Filter by gene if specified
    if gene:
        result = filter_eqtl_by_gene(result, gene)

    # Filter by region if specified
    if chrom is not None and start is not None and end is not None:
        result = filter_eqtl_by_region(result, chrom, start, end, pos_col=pos_col)

    # Add -log10(p) column
    result["neglog10p"] = -np.log10(result[p_col].clip(lower=1e-300))

    return result


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


def calculate_colocalization_overlap(
    gwas_df: pd.DataFrame,
    eqtl_df: pd.DataFrame,
    gwas_pos_col: str = "ps",
    eqtl_pos_col: str = "pos",
    gwas_p_col: str = "p_wald",
    eqtl_p_col: str = "p_value",
    p_threshold: float = 1e-5,
) -> pd.DataFrame:
    """Find SNPs significant in both GWAS and eQTL.

    Simple overlap analysis - for formal colocalization,
    use dedicated tools like coloc or eCAVIAR.

    Args:
        gwas_df: GWAS results DataFrame.
        eqtl_df: eQTL results DataFrame.
        gwas_pos_col: Position column in GWAS data.
        eqtl_pos_col: Position column in eQTL data.
        gwas_p_col: P-value column in GWAS data.
        eqtl_p_col: P-value column in eQTL data.
        p_threshold: P-value threshold for significance.

    Returns:
        DataFrame with overlapping significant SNPs from both datasets.
    """
    # Filter to significant SNPs
    sig_gwas = gwas_df[gwas_df[gwas_p_col] < p_threshold][[gwas_pos_col, gwas_p_col]]
    sig_eqtl = eqtl_df[eqtl_df[eqtl_p_col] < p_threshold][[eqtl_pos_col, eqtl_p_col]]

    # Merge on position
    overlap = sig_gwas.merge(
        sig_eqtl,
        left_on=gwas_pos_col,
        right_on=eqtl_pos_col,
        how="inner",
        suffixes=("_gwas", "_eqtl"),
    )

    logger.info(
        f"Found {len(overlap)} SNPs significant in both GWAS and eQTL "
        f"(p < {p_threshold})"
    )

    return overlap

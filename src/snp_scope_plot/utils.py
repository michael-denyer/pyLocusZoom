"""Utility functions for snp-scope-plot.

Shared helpers used across multiple modules.
"""

from pathlib import Path
from typing import List, Optional, Union

import pandas as pd


class ValidationError(ValueError):
    """Raised when input validation fails."""


def normalize_chrom(chrom: Union[int, str]) -> str:
    """Normalize chromosome identifier by removing 'chr' prefix.

    Args:
        chrom: Chromosome as integer (1, 2, ...) or string ("chr1", "1").

    Returns:
        String without 'chr' prefix (e.g., "1", "X").

    Example:
        >>> normalize_chrom(1)
        '1'
        >>> normalize_chrom("chr1")
        '1'
        >>> normalize_chrom("chrX")
        'X'
    """
    return str(chrom).replace("chr", "")


def validate_dataframe(
    df: pd.DataFrame,
    required_cols: List[str],
    name: str = "DataFrame",
) -> None:
    """Validate that a DataFrame has required columns.

    Args:
        df: DataFrame to validate.
        required_cols: List of required column names.
        name: Name for error messages (e.g., "gwas_df", "genes_df").

    Raises:
        ValidationError: If required columns are missing.

    Example:
        >>> validate_dataframe(df, ["chr", "start", "end"], "genes_df")
    """
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        available = list(df.columns)
        raise ValidationError(
            f"{name} missing required columns: {missing}. "
            f"Available columns: {available}"
        )


def validate_gwas_df(
    df: pd.DataFrame,
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: Optional[str] = None,
) -> None:
    """Validate GWAS results DataFrame.

    Args:
        df: GWAS results DataFrame.
        pos_col: Column name for position.
        p_col: Column name for p-values.
        rs_col: Column name for SNP IDs (optional).

    Raises:
        ValidationError: If required columns are missing.
    """
    required = [pos_col, p_col]
    if rs_col:
        required.append(rs_col)
    validate_dataframe(df, required, "gwas_df")


def validate_genes_df(df: pd.DataFrame) -> None:
    """Validate gene annotations DataFrame.

    Args:
        df: Gene annotations DataFrame.

    Raises:
        ValidationError: If required columns are missing.
    """
    validate_dataframe(df, ["chr", "start", "end", "gene_name"], "genes_df")


def validate_plink_files(bfile_path: Union[str, Path]) -> Path:
    """Validate that PLINK binary fileset exists.

    Checks for .bed, .bim, and .fam files.

    Args:
        bfile_path: Path prefix for PLINK files (without extension).

    Returns:
        Path object if files exist.

    Raises:
        ValidationError: If any PLINK files are missing.
    """
    path = Path(bfile_path)
    missing = []
    for ext in [".bed", ".bim", ".fam"]:
        if not path.with_suffix(ext).exists():
            missing.append(ext)

    if missing:
        raise ValidationError(
            f"PLINK files missing for {path}: {missing}. "
            f"Expected: {path}.bed, {path}.bim, {path}.fam"
        )
    return path

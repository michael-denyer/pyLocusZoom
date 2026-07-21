"""Validation schemas for loaded data.

Provides validation for GWAS, eQTL, fine-mapping, and gene annotation
DataFrames to ensure data quality before plotting.

These are the strict, load-time schemas. Plot-time validation in
``utils.validate_*`` is deliberately more permissive; see the two-tier note in
``CONTEXT.md``.
"""

from pathlib import Path
from typing import Union

import pandas as pd

from .exceptions import LoaderValidationError
from .validation import DataFrameValidator


def _loader_validator(df: pd.DataFrame, name: str) -> DataFrameValidator:
    """Start a validation chain that reports as a loader failure."""
    return DataFrameValidator(df, name, error_class=LoaderValidationError)


# =============================================================================
# GWAS Validation
# =============================================================================


def validate_gwas_dataframe(
    df: pd.DataFrame,
    pos_col: str = "ps",
    p_col: str = "p_wald",
) -> pd.DataFrame:
    """Validate a GWAS DataFrame.

    Args:
        df: DataFrame to validate.
        pos_col: Column name for position.
        p_col: Column name for p-value.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    (
        _loader_validator(df, "GWAS")
        .require_columns([pos_col, p_col])
        .require_numeric([pos_col, p_col])
        .require_not_null([pos_col, p_col])
        .require_range(pos_col, min_val=0, exclusive_min=True)
        .require_range(p_col, min_val=0, max_val=1, exclusive_min=True)
        .validate()
    )
    return df


# =============================================================================
# eQTL Validation
# =============================================================================


def validate_eqtl_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Validate an eQTL DataFrame.

    Args:
        df: DataFrame to validate.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    (
        _loader_validator(df, "eQTL")
        .require_columns(["pos", "p_value", "gene"])
        .require_numeric(["pos", "p_value"])
        .require_range("pos", min_val=0, exclusive_min=True)
        .require_range("p_value", min_val=0, max_val=1, exclusive_min=True)
        .validate()
    )
    return df


# =============================================================================
# Fine-mapping Validation
# =============================================================================


def validate_finemapping_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Validate a fine-mapping DataFrame.

    Args:
        df: DataFrame to validate.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    (
        _loader_validator(df, "Fine-mapping")
        .require_columns(["pos", "pip"])
        .require_numeric(["pos", "pip"])
        .require_range("pos", min_val=0, exclusive_min=True)
        .require_range("pip", min_val=0, max_val=1)
        .validate()
    )
    return df


# =============================================================================
# Gene Annotation Validation
# =============================================================================


def validate_genes_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Validate a genes DataFrame.

    Args:
        df: DataFrame to validate.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    (
        _loader_validator(df, "Gene annotation")
        .require_columns(["chr", "start", "end", "gene_name"])
        .require_numeric(["start", "end"])
        .require_range("start", min_val=0)
        .require_ordering("start", "end")
        .validate()
    )
    return df


# =============================================================================
# File Path Validation
# =============================================================================


def validate_file_path(filepath: Union[str, Path]) -> Path:
    """Validate that a file path exists and is readable.

    Args:
        filepath: Path to validate.

    Returns:
        Validated Path object.

    Raises:
        LoaderValidationError: If file doesn't exist or isn't readable.
    """
    path = Path(filepath)

    if not path.exists():
        raise LoaderValidationError(f"File not found: {path}")

    if not path.is_file():
        raise LoaderValidationError(f"Not a file: {path}")

    return path

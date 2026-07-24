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
from .validation import ColumnSpec, RangeRule, check

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
    check(
        df,
        ColumnSpec(
            name="GWAS",
            required=(pos_col, p_col),
            numeric=(pos_col, p_col),
            not_null=(pos_col, p_col),
            ranges=(RangeRule(pos_col, min_val=0, exclusive_min=True),),
            pvalue=p_col,
            error_class=LoaderValidationError,
        ),
    )
    return df


# =============================================================================
# eQTL Validation
# =============================================================================

_EQTL_SPEC = ColumnSpec(
    name="eQTL",
    required=("pos", "p_value", "gene"),
    numeric=("pos", "p_value"),
    not_null=("pos", "p_value"),
    ranges=(RangeRule("pos", min_val=0, exclusive_min=True),),
    pvalue="p_value",
    error_class=LoaderValidationError,
)


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
    check(df, _EQTL_SPEC)
    return df


# =============================================================================
# Fine-mapping Validation
# =============================================================================

_FINEMAPPING_SPEC = ColumnSpec(
    name="Fine-mapping",
    required=("pos", "pip"),
    numeric=("pos", "pip"),
    not_null=("pos", "pip"),
    ranges=(
        RangeRule("pos", min_val=0, exclusive_min=True),
        RangeRule("pip", min_val=0, max_val=1),
    ),
    error_class=LoaderValidationError,
)


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
    check(df, _FINEMAPPING_SPEC)
    return df


# =============================================================================
# Gene Annotation Validation
# =============================================================================

_GENES_SPEC = ColumnSpec(
    name="Gene annotation",
    required=("chr", "start", "end", "gene_name"),
    numeric=("start", "end"),
    not_null=("start", "end"),
    ranges=(RangeRule("start", min_val=0),),
    ordering=(("start", "end"),),
    error_class=LoaderValidationError,
)


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
    check(df, _GENES_SPEC)
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

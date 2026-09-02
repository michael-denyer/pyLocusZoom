"""Validation schemas for loaded data.

Provides validation for GWAS, eQTL, fine-mapping, and gene annotation
DataFrames to ensure data quality before plotting.

These are the strict, load-time schemas. The permissive plot-time tier is
``validation.gwas_spec(..., strict=False)`` and ``validation.GENES_SPEC``; see
the two-tier note in ``CONTEXT.md``.
"""

import pandas as pd

from .exceptions import LoaderValidationError
from .validation import ColumnSpec, RangeRule, check, gwas_spec

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
    check(df, gwas_spec(pos_col, p_col, strict=True))
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

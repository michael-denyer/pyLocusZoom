"""Colocalization data validation and preparation.

Validates GWAS and eQTL DataFrames for colocalization scatter plots.
"""

from typing import Optional

import pandas as pd

from .validation import ColumnSpec, check


def _validate_coloc_df(
    df: pd.DataFrame,
    df_name: str,
    pos_col: str,
    p_col: str,
    rs_col: Optional[str] = None,
) -> None:
    """Validate DataFrame for colocalization plot.

    Args:
        df: Results DataFrame.
        df_name: Name for error messages (e.g., "GWAS DataFrame").
        pos_col: Column name for genomic positions.
        p_col: Column name for p-values.
        rs_col: Optional column name for SNP IDs.

    Raises:
        ValidationError: If required columns are missing or have invalid types.
    """
    required_cols = [pos_col, p_col]
    if rs_col is not None:
        required_cols.append(rs_col)

    check(
        df,
        ColumnSpec(
            name=df_name,
            required=tuple(required_cols),
            numeric=(pos_col, p_col),
            not_null=(p_col,),
            pvalue=p_col,
        ),
    )


def validate_coloc_gwas_df(
    df: pd.DataFrame,
    pos_col: str,
    p_col: str,
    rs_col: Optional[str] = None,
) -> None:
    """Validate GWAS DataFrame for colocalization plot.

    Args:
        df: GWAS results DataFrame.
        pos_col: Column name for genomic positions.
        p_col: Column name for p-values.
        rs_col: Optional column name for SNP IDs.

    Raises:
        ValidationError: If required columns are missing or have invalid types.
    """
    _validate_coloc_df(df, "GWAS DataFrame", pos_col, p_col, rs_col)


def validate_coloc_eqtl_df(
    df: pd.DataFrame,
    pos_col: str,
    p_col: str,
    rs_col: Optional[str] = None,
) -> None:
    """Validate eQTL DataFrame for colocalization plot.

    Args:
        df: eQTL results DataFrame.
        pos_col: Column name for genomic positions.
        p_col: Column name for p-values.
        rs_col: Optional column name for SNP IDs.

    Raises:
        ValidationError: If required columns are missing or have invalid types.
    """
    _validate_coloc_df(df, "eQTL DataFrame", pos_col, p_col, rs_col)

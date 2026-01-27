"""Fine-mapping/SuSiE data handling for pyLocusZoom.

Provides utilities for loading, validating, and preparing statistical
fine-mapping results (SuSiE, FINEMAP, etc.) for visualization.
"""

from typing import List, Optional

import pandas as pd

from .logging import logger
from .utils import ValidationError
from .validation import DataFrameValidator

# Required columns for fine-mapping data
REQUIRED_FINEMAPPING_COLS = ["pos", "pip"]
OPTIONAL_FINEMAPPING_COLS = ["rs", "cs", "cs_id", "effect", "se"]


class FinemappingValidationError(ValueError):
    """Raised when fine-mapping DataFrame validation fails."""

    pass


def validate_finemapping_df(
    df: pd.DataFrame,
    pos_col: str = "pos",
    pip_col: str = "pip",
) -> None:
    """Validate fine-mapping DataFrame has required columns.

    Args:
        df: Fine-mapping DataFrame to validate.
        pos_col: Column name for genomic position.
        pip_col: Column name for posterior inclusion probability.

    Raises:
        FinemappingValidationError: If required columns are missing.
    """
    try:
        (
            DataFrameValidator(df, "Fine-mapping DataFrame")
            .require_columns([pos_col, pip_col])
            .require_numeric([pip_col])
            .require_range(pip_col, min_val=0, max_val=1)
            .validate()
        )
    except ValidationError as e:
        raise FinemappingValidationError(str(e)) from e


def filter_finemapping_by_region(
    df: pd.DataFrame,
    chrom: int,
    start: int,
    end: int,
    pos_col: str = "pos",
    chrom_col: Optional[str] = "chr",
) -> pd.DataFrame:
    """Filter fine-mapping data to a genomic region.

    Args:
        df: Fine-mapping DataFrame.
        chrom: Chromosome number.
        start: Start position.
        end: End position.
        pos_col: Column name for position.
        chrom_col: Column name for chromosome (if present).

    Returns:
        Filtered DataFrame containing only variants in the region.
    """
    mask = (df[pos_col] >= start) & (df[pos_col] <= end)

    # Filter by chromosome if column exists
    if chrom_col and chrom_col in df.columns:
        chrom_str = str(chrom).replace("chr", "")
        df_chrom = df[chrom_col].astype(str).str.replace("chr", "", regex=False)
        mask = mask & (df_chrom == chrom_str)

    filtered = df[mask].copy()
    logger.debug(
        f"Filtered fine-mapping data to {len(filtered)} variants in region "
        f"chr{chrom}:{start}-{end}"
    )
    return filtered


def get_credible_sets(
    df: pd.DataFrame,
    cs_col: str = "cs",
) -> List[int]:
    """Get list of unique credible set IDs.

    Args:
        df: Fine-mapping DataFrame.
        cs_col: Column containing credible set assignments.

    Returns:
        Sorted list of unique credible set IDs (excluding 0/NA).
    """
    if cs_col not in df.columns:
        return []
    # Filter out variants not in a credible set (typically cs=0 or NA)
    cs_values = df[cs_col].dropna()
    cs_values = cs_values[cs_values != 0]
    return sorted(cs_values.unique().tolist())


def filter_by_credible_set(
    df: pd.DataFrame,
    cs_id: int,
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Filter to variants in a specific credible set.

    Args:
        df: Fine-mapping DataFrame.
        cs_id: Credible set ID to filter for.
        cs_col: Column containing credible set assignments.

    Returns:
        Filtered DataFrame containing only variants in the credible set.
    """
    if cs_col not in df.columns:
        raise FinemappingValidationError(
            f"Cannot filter by credible set: column '{cs_col}' not found. "
            f"Available columns: {list(df.columns)}"
        )
    return df[df[cs_col] == cs_id].copy()


def prepare_finemapping_for_plotting(
    df: pd.DataFrame,
    pos_col: str = "pos",
    pip_col: str = "pip",
    chrom: Optional[int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
) -> pd.DataFrame:
    """Prepare fine-mapping data for plotting.

    Validates, filters, and sorts data for plotting as a line or scatter.

    Args:
        df: Raw fine-mapping DataFrame.
        pos_col: Column name for position.
        pip_col: Column name for PIP.
        chrom: Optional chromosome for region filtering.
        start: Optional start position for region filtering.
        end: Optional end position for region filtering.

    Returns:
        Prepared DataFrame sorted by position.
    """
    validate_finemapping_df(df, pos_col=pos_col, pip_col=pip_col)

    result = df.copy()

    # Filter by region if specified
    if chrom is not None and start is not None and end is not None:
        result = filter_finemapping_by_region(
            result, chrom, start, end, pos_col=pos_col
        )

    # Sort by position for line plotting
    result = result.sort_values(pos_col)

    return result


def get_top_pip_variants(
    df: pd.DataFrame,
    n: int = 5,
    pip_col: str = "pip",
    pip_threshold: float = 0.0,
) -> pd.DataFrame:
    """Get top variants by posterior inclusion probability.

    Args:
        df: Fine-mapping DataFrame.
        n: Number of top variants to return.
        pip_col: Column containing PIP values.
        pip_threshold: Minimum PIP threshold.

    Returns:
        DataFrame with top N variants by PIP.
    """
    filtered = df[df[pip_col] >= pip_threshold]
    return filtered.nlargest(n, pip_col)


def calculate_credible_set_coverage(
    df: pd.DataFrame,
    cs_col: str = "cs",
    pip_col: str = "pip",
) -> dict:
    """Calculate cumulative PIP for each credible set.

    Args:
        df: Fine-mapping DataFrame.
        cs_col: Column containing credible set assignments.
        pip_col: Column containing PIP values.

    Returns:
        Dictionary mapping credible set ID to cumulative PIP.
    """
    if cs_col not in df.columns:
        return {}

    coverage = {}
    for cs_id in get_credible_sets(df, cs_col):
        cs_data = filter_by_credible_set(df, cs_id, cs_col)
        coverage[cs_id] = cs_data[pip_col].sum()

    return coverage

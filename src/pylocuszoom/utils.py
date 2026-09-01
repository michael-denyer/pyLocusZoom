"""Utility functions for pyLocusZoom.

Shared helpers used across multiple modules.

"""

import os
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, List, Optional, Union

import pandas as pd

from .exceptions import ValidationError

if TYPE_CHECKING:
    from pyspark.sql import DataFrame as SparkDataFrame

# Type alias for DataFrames (pandas or PySpark)
DataFrameLike = Union[pd.DataFrame, "SparkDataFrame"]


def _platform_cache_base() -> Path:
    """Return the platform-appropriate pyLocusZoom cache root.

    Single source of truth for cache placement. Recombination maps and Ensembl
    annotations both live under this root; each appends its own leaf
    (``recombination_maps`` / ``ensembl``).

    Resolution order:
    - Windows: %LOCALAPPDATA%/pylocuszoom (or ~/AppData/Local/pylocuszoom if unset)
    - Databricks (``/dbfs`` present): /dbfs/FileStore/reference_data
    - macOS/Linux: $XDG_CACHE_HOME/pylocuszoom (or ~/.cache/pylocuszoom)
    """
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif os.path.exists("/dbfs"):
        return Path("/dbfs/FileStore/reference_data")
    else:
        xdg_cache = os.environ.get("XDG_CACHE_HOME")
        base = Path(xdg_cache) if xdg_cache else Path.home() / ".cache"

    return base / "pylocuszoom"


def is_spark_dataframe(df: Any) -> bool:
    """Check if object is a PySpark DataFrame.

    Args:
        df: Object to check.

    Returns:
        True if PySpark DataFrame, False otherwise.
    """
    # Check class name to avoid importing pyspark
    return type(df).__name__ == "DataFrame" and type(df).__module__.startswith(
        "pyspark"
    )


def to_pandas(
    df: DataFrameLike,
    sample_size: Optional[int] = None,
) -> pd.DataFrame:
    """Convert DataFrame-like object to pandas DataFrame.

    Supports pandas DataFrames (returned as-is) and PySpark DataFrames
    (converted to pandas). For large PySpark DataFrames, use sample_size
    to limit the data transferred.

    Args:
        df: pandas DataFrame or PySpark DataFrame.
        sample_size: For PySpark, limit to this many rows. If None,
            converts entire DataFrame (may be slow for large data).

    Returns:
        pandas DataFrame.

    Raises:
        TypeError: If df is not a supported DataFrame type.

    Example:
        >>> # PySpark DataFrame
        >>> pdf = to_pandas(spark_df, sample_size=100000)
        >>>
        >>> # pandas DataFrame (passthrough)
        >>> pdf = to_pandas(pandas_df)
    """
    if isinstance(df, pd.DataFrame):
        return df

    if is_spark_dataframe(df):
        if sample_size is not None:
            # Sample to limit data transfer
            total = df.count()
            if total > sample_size:
                fraction = sample_size / total
                df = df.sample(fraction=fraction, seed=42)
        return df.toPandas()

    # Try pandas conversion as fallback
    if hasattr(df, "to_pandas"):
        return df.to_pandas()
    if hasattr(df, "toPandas"):
        return df.toPandas()

    raise TypeError(
        f"Unsupported DataFrame type: {type(df).__name__}. "
        f"Expected pandas.DataFrame or pyspark.sql.DataFrame"
    )


# Different spellings of the same assembly, keyed by assembly_token output.
ASSEMBLY_SYNONYMS: dict[str, str] = {
    # Canine
    "canfam31": "canfam3",
    "canfam40": "canfam4",
    "uucfamgsd10": "canfam4",
    "roscfam10": "roscfam1",
    # Feline
    "felcat90": "felcat9",
    "feliscatus90": "felcat9",
    "fcatusfca126mat10": "fca126",
    # Human
    "hg19": "grch37",
    "grch37p13": "grch37",
    "hg38": "grch38",
    "grch38p14": "grch38",
    # Mouse
    "mm39": "grcm39",
}


def assembly_token(name: str) -> str:
    """Reduce an assembly or build name to a comparable token.

    Strips punctuation and case, then folds known synonyms together, so
    ``"CanFam4.0"``, ``"canfam4"`` and ``"UU_Cfam_GSD_1.0"`` all compare equal.
    """
    token = re.sub(r"[^a-z0-9]", "", name.lower())
    return ASSEMBLY_SYNONYMS.get(token, token)


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


def filter_by_region(
    df: pd.DataFrame,
    region: tuple,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
) -> pd.DataFrame:
    """Filter DataFrame to genomic region with inclusive bounds.

    Filters rows where position is within [start, end] (inclusive).
    If chrom_col exists in DataFrame, also filters by chromosome.
    Chromosome comparison normalizes types (int/str, chr prefix).

    Args:
        df: DataFrame to filter.
        region: Tuple of (chrom, start, end) defining the region.
        chrom_col: Column name for chromosome (default: "chrom").
            If column doesn't exist, filters by position only.
        pos_col: Column name for position (default: "pos").

    Returns:
        Filtered DataFrame (copy, not view).

    Raises:
        KeyError: If pos_col is not found in DataFrame.

    Example:
        >>> filtered = filter_by_region(df, region=(1, 1000000, 2000000))
        >>> filtered = filter_by_region(df, region=("chr1", 1e6, 2e6), pos_col="position")
    """
    chrom, start, end = region

    # Validate position column exists
    if pos_col not in df.columns:
        raise KeyError(
            f"Position column '{pos_col}' not found in DataFrame. "
            f"Available columns: {list(df.columns)}"
        )

    # Position filtering (inclusive bounds)
    mask = (df[pos_col] >= start) & (df[pos_col] <= end)

    # Chromosome filtering (if column exists)
    if chrom_col in df.columns:
        chrom_normalized = normalize_chrom(chrom)
        df_chrom_normalized = (
            df[chrom_col].astype(str).str.replace("chr", "", regex=False)
        )
        mask = mask & (df_chrom_normalized == chrom_normalized)

    return df[mask].copy()


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

    This is the permissive plot-time tier: columns and non-emptiness only. The
    p-value policy is deferred to ``_data.prepare_pvalue_data``, which keeps an
    exact zero under the Manhattan convention. The strict ``(0, 1]`` tier lives
    in ``schemas.validate_gwas_dataframe`` and runs at load time. The split is
    deliberate: tightening this function would reject Manhattan input that
    plots correctly today.

    Args:
        df: GWAS results DataFrame.
        pos_col: Column name for position.
        p_col: Column name for p-values.
        rs_col: Column name for SNP IDs (optional).

    Raises:
        ValidationError: If required columns are missing or DataFrame is empty.
    """
    if df.empty:
        raise ValidationError("gwas_df is empty — no rows to plot")

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

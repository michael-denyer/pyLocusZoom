"""Utility functions for pyLocusZoom.

Shared helpers used across multiple modules.

"""

import os
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, Union

import pandas as pd

from .exceptions import ValidationError
from .schemas import Canonical

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

    Every public plot method calls this on the frames it is given, so the
    library's own boundary is the only place a Spark frame is collected.

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
        ValidationError: If df is not a supported DataFrame type.

    Example:
        >>> # PySpark DataFrame
        >>> pdf = to_pandas(spark_df, sample_size=100000)
        >>>
        >>> # pandas DataFrame (passthrough)
        >>> pdf = to_pandas(pandas_df)
    """
    if isinstance(df, pd.DataFrame):
        return df

    if is_spark_dataframe(df) and sample_size is not None:
        # Sample to limit data transfer
        total = df.count()
        if total > sample_size:
            fraction = sample_size / total
            df = df.sample(fraction=fraction, seed=42)

    # toPandas() is the Spark contract, and the only one this function needs:
    # Databricks Connect and other wrappers satisfy it without being importable
    # as pyspark, so is_spark_dataframe alone would turn them away.
    if hasattr(df, "toPandas"):
        return df.toPandas()

    raise ValidationError(
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


# A leading "chr" in any case; one later in the name (chrUn_chr5) is kept.
_CHR_PREFIX = re.compile(r"^chr", re.IGNORECASE)


def normalize_chrom(chrom: Union[int, str]) -> str:
    """Normalize chromosome identifier by removing a leading 'chr', in any case.

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
        >>> normalize_chrom("Chr1")
        '1'
    """
    return _CHR_PREFIX.sub("", str(chrom))


def normalize_chrom_series(chroms: pd.Series) -> pd.Series:
    """Normalize a column of chromosome identifiers by removing a 'chr' prefix.

    The frame-level companion to :func:`normalize_chrom`, so a column of mixed
    integer and ``"chr1"`` spellings compares equal to a normalized scalar.

    Args:
        chroms: Column of chromosome identifiers, of any dtype.

    Returns:
        A string Series without 'chr' prefixes.

    Example:
        >>> normalize_chrom_series(pd.Series([1, "chr2", "chrX"])).tolist()
        ['1', '2', 'X']
        >>> normalize_chrom_series(pd.Series([1.0, 2.0])).tolist()
        ['1', '2']
    """
    text = chroms.astype(str)
    if pd.api.types.is_float_dtype(chroms):
        # pandas reads an integer column holding a NaN as float; 1.0 means "1".
        integral = chroms.notna() & (chroms % 1 == 0)
        as_int = chroms.where(integral, 0).astype("int64").astype(str)
        text = text.mask(integral, as_int)
    return text.str.replace(_CHR_PREFIX, "", regex=True)


def filter_by_region(
    df: pd.DataFrame,
    region: tuple,
    chrom_col: str | None = Canonical.CHROM,
    pos_col: str = Canonical.POS,
) -> pd.DataFrame:
    """Filter DataFrame to genomic region with inclusive bounds.

    Filters rows where position is within [start, end] (inclusive) and the
    chromosome matches, compared through :func:`normalize_chrom_series`.

    Args:
        df: DataFrame to filter.
        region: Tuple of (chrom, start, end) defining the region.
        chrom_col: Column name for chromosome (default: "chr"). None selects
            by position only, for a frame already scoped to one chromosome.
        pos_col: Column name for position (default: "pos").

    Returns:
        Filtered DataFrame (copy, not view).

    Raises:
        ValidationError: If pos_col, or chrom_col when not None, is not found
            in DataFrame.

    Example:
        >>> filtered = filter_by_region(df, region=(1, 1000000, 2000000))
        >>> filtered = filter_by_region(df, region=("chr1", 1e6, 2e6), pos_col="position")
    """
    chrom, start, end = region

    # Validate position column exists
    if pos_col not in df.columns:
        raise ValidationError(
            f"Position column '{pos_col}' not found in DataFrame. "
            f"Available columns: {list(df.columns)}"
        )

    # Position filtering (inclusive bounds)
    mask = (df[pos_col] >= start) & (df[pos_col] <= end)

    if chrom_col is not None:
        if chrom_col not in df.columns:
            raise ValidationError(
                f"Chromosome column '{chrom_col}' not found in DataFrame. "
                f"Available columns: {list(df.columns)}. Pass chrom_col=None "
                "to select by position only."
            )
        mask = mask & (normalize_chrom_series(df[chrom_col]) == normalize_chrom(chrom))

    return df[mask].copy()

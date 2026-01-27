"""Pydantic validation schemas for loaded data.

Provides validation models for GWAS, eQTL, fine-mapping, and gene annotation
DataFrames to ensure data quality before plotting.
"""

from pathlib import Path
from typing import Optional, Union

import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator, model_validator


class LoaderValidationError(Exception):
    """Raised when loaded data fails validation."""

    pass


# =============================================================================
# GWAS Validation
# =============================================================================


class GWASRowModel(BaseModel):
    """Validation model for a single GWAS row."""

    model_config = ConfigDict(extra="allow")

    ps: int
    p_wald: float
    rs: Optional[str] = None
    chr: Optional[Union[str, int]] = None

    @field_validator("ps")
    @classmethod
    def position_positive(cls, v: int) -> int:
        """Position must be positive."""
        if v <= 0:
            raise ValueError(f"Position must be positive, got {v}")
        return v

    @field_validator("p_wald")
    @classmethod
    def pvalue_in_range(cls, v: float) -> float:
        """P-value must be between 0 and 1."""
        if not (0 < v <= 1):
            raise ValueError(f"P-value must be in range (0, 1], got {v}")
        return v


def validate_gwas_dataframe(
    df: pd.DataFrame,
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
    strict: bool = False,
) -> pd.DataFrame:
    """Validate a GWAS DataFrame.

    Args:
        df: DataFrame to validate.
        pos_col: Column name for position.
        p_col: Column name for p-value.
        rs_col: Column name for SNP ID.
        strict: If True, validate every row. If False (default), validate schema only.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    errors = []

    # Check required columns exist
    if pos_col not in df.columns:
        errors.append(f"Missing required column: '{pos_col}'")
    if p_col not in df.columns:
        errors.append(f"Missing required column: '{p_col}'")

    if errors:
        raise LoaderValidationError(
            "GWAS validation failed:\n  - " + "\n  - ".join(errors)
        )

    # Check data types (must be numeric for range checks)
    pos_is_numeric = pd.api.types.is_numeric_dtype(df[pos_col])
    p_is_numeric = pd.api.types.is_numeric_dtype(df[p_col])

    if not pos_is_numeric:
        errors.append(f"Column '{pos_col}' must be numeric, got {df[pos_col].dtype}")

    if not p_is_numeric:
        errors.append(f"Column '{p_col}' must be numeric, got {df[p_col].dtype}")

    # Only check value ranges if columns are numeric (avoid confusing errors)
    if pos_is_numeric:
        if (df[pos_col] <= 0).any():
            n_invalid = (df[pos_col] <= 0).sum()
            errors.append(f"Column '{pos_col}' has {n_invalid} non-positive values")

        if df[pos_col].isna().any():
            n_na = df[pos_col].isna().sum()
            errors.append(f"Column '{pos_col}' has {n_na} missing values")

    if p_is_numeric:
        if ((df[p_col] <= 0) | (df[p_col] > 1)).any():
            n_invalid = ((df[p_col] <= 0) | (df[p_col] > 1)).sum()
            errors.append(
                f"Column '{p_col}' has {n_invalid} values outside range (0, 1]"
            )

        if df[p_col].isna().any():
            n_na = df[p_col].isna().sum()
            errors.append(f"Column '{p_col}' has {n_na} missing values")

    if errors:
        raise LoaderValidationError(
            "GWAS validation failed:\n  - " + "\n  - ".join(errors)
        )

    return df


# =============================================================================
# eQTL Validation
# =============================================================================


class EQTLRowModel(BaseModel):
    """Validation model for a single eQTL row."""

    model_config = ConfigDict(extra="allow")

    pos: int
    p_value: float
    gene: str
    effect: Optional[float] = None

    @field_validator("pos")
    @classmethod
    def position_positive(cls, v: int) -> int:
        """Position must be positive."""
        if v <= 0:
            raise ValueError(f"Position must be positive, got {v}")
        return v

    @field_validator("p_value")
    @classmethod
    def pvalue_in_range(cls, v: float) -> float:
        """P-value must be between 0 and 1."""
        if not (0 < v <= 1):
            raise ValueError(f"P-value must be in range (0, 1], got {v}")
        return v


def validate_eqtl_dataframe(
    df: pd.DataFrame,
    strict: bool = False,
) -> pd.DataFrame:
    """Validate an eQTL DataFrame.

    Args:
        df: DataFrame to validate.
        strict: If True, validate every row.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    errors = []

    # Check required columns
    required = ["pos", "p_value", "gene"]
    for col in required:
        if col not in df.columns:
            errors.append(f"Missing required column: '{col}'")

    if errors:
        raise LoaderValidationError(
            "eQTL validation failed:\n  - " + "\n  - ".join(errors)
        )

    # Check data types and ranges
    if not pd.api.types.is_numeric_dtype(df["pos"]):
        errors.append(f"Column 'pos' must be numeric, got {df['pos'].dtype}")
    elif (df["pos"] <= 0).any():
        n_invalid = (df["pos"] <= 0).sum()
        errors.append(f"Column 'pos' has {n_invalid} non-positive values")

    if not pd.api.types.is_numeric_dtype(df["p_value"]):
        errors.append(f"Column 'p_value' must be numeric, got {df['p_value'].dtype}")
    elif ((df["p_value"] <= 0) | (df["p_value"] > 1)).any():
        n_invalid = ((df["p_value"] <= 0) | (df["p_value"] > 1)).sum()
        errors.append(f"Column 'p_value' has {n_invalid} values outside range (0, 1]")

    if errors:
        raise LoaderValidationError(
            "eQTL validation failed:\n  - " + "\n  - ".join(errors)
        )

    return df


# =============================================================================
# Fine-mapping Validation
# =============================================================================


class FinemappingRowModel(BaseModel):
    """Validation model for a single fine-mapping row."""

    model_config = ConfigDict(extra="allow")

    pos: int
    pip: float
    cs: Optional[int] = None

    @field_validator("pos")
    @classmethod
    def position_positive(cls, v: int) -> int:
        """Position must be positive."""
        if v <= 0:
            raise ValueError(f"Position must be positive, got {v}")
        return v

    @field_validator("pip")
    @classmethod
    def pip_in_range(cls, v: float) -> float:
        """PIP must be between 0 and 1."""
        if not (0 <= v <= 1):
            raise ValueError(f"PIP must be in range [0, 1], got {v}")
        return v


def validate_finemapping_dataframe(
    df: pd.DataFrame,
    cs_col: str = "cs",
    strict: bool = False,
) -> pd.DataFrame:
    """Validate a fine-mapping DataFrame.

    Args:
        df: DataFrame to validate.
        cs_col: Column name for credible set.
        strict: If True, validate every row.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    errors = []

    # Check required columns
    if "pos" not in df.columns:
        errors.append("Missing required column: 'pos'")
    if "pip" not in df.columns:
        errors.append("Missing required column: 'pip'")

    if errors:
        raise LoaderValidationError(
            "Fine-mapping validation failed:\n  - " + "\n  - ".join(errors)
        )

    # Check data types and ranges
    if not pd.api.types.is_numeric_dtype(df["pos"]):
        errors.append(f"Column 'pos' must be numeric, got {df['pos'].dtype}")
    elif (df["pos"] <= 0).any():
        n_invalid = (df["pos"] <= 0).sum()
        errors.append(f"Column 'pos' has {n_invalid} non-positive values")

    if not pd.api.types.is_numeric_dtype(df["pip"]):
        errors.append(f"Column 'pip' must be numeric, got {df['pip'].dtype}")
    elif ((df["pip"] < 0) | (df["pip"] > 1)).any():
        n_invalid = ((df["pip"] < 0) | (df["pip"] > 1)).sum()
        errors.append(f"Column 'pip' has {n_invalid} values outside range [0, 1]")

    if errors:
        raise LoaderValidationError(
            "Fine-mapping validation failed:\n  - " + "\n  - ".join(errors)
        )

    return df


# =============================================================================
# Gene Annotation Validation
# =============================================================================


class GeneRowModel(BaseModel):
    """Validation model for a single gene annotation row."""

    model_config = ConfigDict(extra="allow")

    chr: Union[str, int]
    start: int
    end: int
    gene_name: str
    strand: Optional[str] = None

    @field_validator("start", "end")
    @classmethod
    def position_positive(cls, v: int) -> int:
        """Position must be positive."""
        if v < 0:
            raise ValueError(f"Position must be non-negative, got {v}")
        return v

    @model_validator(mode="after")
    def start_before_end(self):
        """Start must be <= end."""
        if self.start > self.end:
            raise ValueError(f"Start ({self.start}) must be <= end ({self.end})")
        return self


def validate_genes_dataframe(
    df: pd.DataFrame,
    strict: bool = False,
) -> pd.DataFrame:
    """Validate a genes DataFrame.

    Args:
        df: DataFrame to validate.
        strict: If True, validate every row.

    Returns:
        Validated DataFrame.

    Raises:
        LoaderValidationError: If validation fails.
    """
    errors = []

    # Check required columns
    required = ["chr", "start", "end", "gene_name"]
    for col in required:
        if col not in df.columns:
            errors.append(f"Missing required column: '{col}'")

    if errors:
        raise LoaderValidationError(
            "Gene annotation validation failed:\n  - " + "\n  - ".join(errors)
        )

    # Check data types
    start_is_numeric = pd.api.types.is_numeric_dtype(df["start"])
    end_is_numeric = pd.api.types.is_numeric_dtype(df["end"])

    if not start_is_numeric:
        errors.append(f"Column 'start' must be numeric, got {df['start'].dtype}")

    if not end_is_numeric:
        errors.append(f"Column 'end' must be numeric, got {df['end'].dtype}")

    # Only check ranges if columns are numeric (avoid confusing errors)
    if start_is_numeric:
        if (df["start"] < 0).any():
            n_invalid = (df["start"] < 0).sum()
            errors.append(f"Column 'start' has {n_invalid} negative values")

    if start_is_numeric and end_is_numeric:
        if (df["end"] < df["start"]).any():
            n_invalid = (df["end"] < df["start"]).sum()
            errors.append(f"Found {n_invalid} genes where end < start")

    if errors:
        raise LoaderValidationError(
            "Gene annotation validation failed:\n  - " + "\n  - ".join(errors)
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

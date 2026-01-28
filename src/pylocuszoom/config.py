"""Pydantic configuration classes for pyLocusZoom plot methods.

This module provides typed, validated configuration objects that replace
the parameter explosion in plot methods. Each config is immutable (frozen)
to prevent accidental modification.

Example:
    >>> from pylocuszoom.config import RegionConfig, DisplayConfig
    >>> region = RegionConfig(chrom=1, start=1000000, end=2000000)
    >>> display = DisplayConfig(snp_labels=False, label_top_n=3)
"""

from typing import Optional, Tuple

from pydantic import BaseModel, ConfigDict, Field, model_validator


class RegionConfig(BaseModel):
    """Genomic region specification.

    Attributes:
        chrom: Chromosome number (must be >= 1).
        start: Start position in base pairs (must be >= 0).
        end: End position in base pairs (must be > start).
    """

    model_config = ConfigDict(frozen=True)

    chrom: int = Field(..., ge=1, description="Chromosome number")
    start: int = Field(..., ge=0, description="Start position (bp)")
    end: int = Field(..., gt=0, description="End position (bp)")

    @model_validator(mode="after")
    def validate_region(self) -> "RegionConfig":
        """Validate that start < end."""
        if self.start >= self.end:
            raise ValueError(f"start ({self.start}) must be < end ({self.end})")
        return self


class ColumnConfig(BaseModel):
    """DataFrame column name mappings for GWAS data.

    Attributes:
        pos_col: Column name for genomic position.
        p_col: Column name for p-value.
        rs_col: Column name for SNP identifier.
    """

    model_config = ConfigDict(frozen=True)

    pos_col: str = Field(default="ps", description="Position column name")
    p_col: str = Field(default="p_wald", description="P-value column name")
    rs_col: str = Field(default="rs", description="SNP ID column name")


class DisplayConfig(BaseModel):
    """Display and visual options for plots.

    Attributes:
        snp_labels: Whether to show SNP labels on plot.
        label_top_n: Number of top SNPs to label.
        show_recombination: Whether to show recombination rate overlay.
        figsize: Figure size as (width, height) in inches.
    """

    model_config = ConfigDict(frozen=True)

    snp_labels: bool = Field(default=True, description="Show SNP labels")
    label_top_n: int = Field(default=5, ge=0, description="Number of top SNPs to label")
    show_recombination: bool = Field(
        default=True, description="Show recombination overlay"
    )
    figsize: Tuple[float, float] = Field(
        default=(12.0, 8.0), description="Figure size (width, height)"
    )


class LDConfig(BaseModel):
    """Linkage disequilibrium configuration.

    Supports three modes:
    1. No LD coloring: All fields None (default)
    2. Pre-computed LD: Provide ld_col for column with R^2 values
    3. Calculate LD: Provide lead_pos and ld_reference_file

    Attributes:
        lead_pos: Position of lead/index SNP to highlight.
        ld_reference_file: Path to PLINK binary fileset for LD calculation.
        ld_col: Column name for pre-computed LD (R^2) values.
    """

    model_config = ConfigDict(frozen=True)

    lead_pos: Optional[int] = Field(default=None, ge=1, description="Lead SNP position")
    ld_reference_file: Optional[str] = Field(
        default=None, description="PLINK binary fileset path"
    )
    ld_col: Optional[str] = Field(
        default=None, description="Pre-computed LD column name"
    )

    @model_validator(mode="after")
    def validate_ld_config(self) -> "LDConfig":
        """Validate LD configuration consistency.

        When ld_reference_file is provided, lead_pos is required to identify
        the index SNP for LD calculation.
        """
        if self.ld_reference_file is not None and self.lead_pos is None:
            raise ValueError("lead_pos is required when ld_reference_file is provided")
        return self


__all__ = ["RegionConfig", "ColumnConfig", "DisplayConfig", "LDConfig"]

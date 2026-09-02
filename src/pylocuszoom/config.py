"""Pydantic configuration classes for pyLocusZoom plot methods.

This module provides typed, validated configuration objects that replace
the parameter explosion in plot methods. Each config is immutable (frozen)
to prevent accidental modification.

Example:
    >>> from pylocuszoom.config import RegionConfig, DisplayConfig, PlotConfig
    >>> region = RegionConfig(chrom=1, start=1000000, end=2000000)
    >>> display = DisplayConfig(snp_labels=False, label_top_n=3)
    >>>
    >>> # Using composite PlotConfig with factory method
    >>> config = PlotConfig.from_kwargs(chrom=1, start=1000000, end=2000000)
"""

from typing import Any, List, Optional, Tuple, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class RegionConfig(BaseModel):
    """Genomic region specification.

    Attributes:
        chrom: Chromosome number or name (e.g., 1, "A1", "X").
        start: Start position in base pairs (1-based, must be >= 1).
        end: End position in base pairs (must be > start).
    """

    model_config = ConfigDict(frozen=True)

    chrom: Union[int, str] = Field(..., description="Chromosome number or name")
    start: int = Field(..., ge=1, description="Start position (bp)")
    end: int = Field(..., gt=0, description="End position (bp)")

    @field_validator("chrom")
    @classmethod
    def validate_chrom(cls, v: Union[int, str]) -> Union[int, str]:
        """Validate chromosome: int must be >= 1, string must be non-empty."""
        if isinstance(v, int) and v < 1:
            raise ValueError(f"Integer chromosome must be >= 1, got {v}")
        if isinstance(v, str) and not v.strip():
            raise ValueError("Chromosome string must not be empty")
        return v

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
        """Validate LD configuration is not contradictory.

        ld_col means LD is pre-computed in a DataFrame column.
        ld_reference_file means LD should be calculated from a PLINK fileset.
        These are mutually exclusive -- using both is ambiguous.
        """
        if self.ld_col is not None and self.ld_reference_file is not None:
            raise ValueError(
                "Cannot specify both ld_col (pre-computed LD) and "
                "ld_reference_file (compute LD). Choose one."
            )
        return self


class PlotConfig(BaseModel):
    """Composite configuration for plot() method.

    Composes all sub-configs into a single validated configuration object.
    Use either direct construction with nested configs, or the from_kwargs()
    factory method for backward compatibility with existing code.

    Attributes:
        region: Genomic region specification (required).
        columns: DataFrame column name mappings.
        display: Display and visual options.
        ld: Linkage disequilibrium configuration.

    Example:
        >>> # Direct construction
        >>> config = PlotConfig(
        ...     region=RegionConfig(chrom=1, start=1000000, end=2000000),
        ...     display=DisplayConfig(snp_labels=False),
        ... )
        >>>
        >>> # Factory method (backward compatible with plot() signature)
        >>> config = PlotConfig.from_kwargs(
        ...     chrom=1, start=1000000, end=2000000,
        ...     snp_labels=False, lead_pos=1500000,
        ... )
    """

    model_config = ConfigDict(frozen=True)

    region: RegionConfig
    columns: ColumnConfig = Field(default_factory=ColumnConfig)
    display: DisplayConfig = Field(default_factory=DisplayConfig)
    ld: LDConfig = Field(default_factory=LDConfig)

    @model_validator(mode="after")
    def validate_ld_requires_lead_pos(self) -> "PlotConfig":
        """Validate that LD reference file has lead_pos for single plots."""
        if self.ld.ld_reference_file is not None and self.ld.lead_pos is None:
            raise ValueError("lead_pos is required when ld_reference_file is provided")
        return self

    @classmethod
    def from_kwargs(cls, **kwargs: Any) -> "PlotConfig":
        """Create a config from the flat keyword arguments a plot method accepts.

        Each keyword is routed to the nested model that declares it; whatever
        remains is a field of ``cls`` itself.

        Raises:
            ValidationError: If parameters are invalid.
        """
        nested = (
            ("region", RegionConfig),
            ("columns", ColumnConfig),
            ("display", DisplayConfig),
            ("ld", LDConfig),
        )
        parts = {
            name: model(
                **{k: kwargs.pop(k) for k in list(kwargs) if k in model.model_fields}
            )
            for name, model in nested
        }
        return cls(**parts, **kwargs)


class StackedPlotConfig(PlotConfig):
    """Composite configuration for plot_stacked() method.

    Extends PlotConfig with list-based parameters for stacked plots.
    Supports multiple lead positions, panel labels, and LD reference files.

    Attributes:
        lead_positions: List of lead SNP positions (one per panel).
        panel_labels: List of panel labels (one per panel).
        ld_reference_files: List of PLINK filesets (one per panel).

    Example:
        >>> config = StackedPlotConfig.from_kwargs(
        ...     chrom=1, start=1000000, end=2000000,
        ...     lead_positions=[1500000, 1600000],
        ...     panel_labels=["Study A", "Study B"],
        ... )
    """

    lead_positions: Optional[List[int]] = Field(
        default=None, description="Lead SNP positions (one per panel)"
    )
    panel_labels: Optional[List[str]] = Field(
        default=None, description="Panel labels (one per panel)"
    )
    ld_reference_files: Optional[List[str]] = Field(
        default=None, description="PLINK filesets (one per panel)"
    )

    @model_validator(mode="after")
    def validate_ld_requires_lead_pos(self) -> "StackedPlotConfig":
        """Validate broadcast LD configuration for stacked plots.

        When ld_reference_file is provided for broadcast, lead_positions must
        be provided to specify the reference SNP for each panel.
        """
        if (
            self.ld.ld_reference_file is not None
            and self.ld.lead_pos is None
            and self.lead_positions is None
        ):
            raise ValueError(
                "lead_positions is required when ld_reference_file is provided "
                "for broadcast (one lead position per panel)"
            )
        return self


class ColocConfig(BaseModel):
    """Configuration for colocalization plot.

    Attributes:
        gwas_p_col: Column name for GWAS p-values.
        eqtl_p_col: Column name for eQTL p-values.
        pos_col: Column name for genomic position.
        rs_col: Optional column name for SNP identifiers.
        ld_col: Optional column name for pre-computed LD values.
        lead_snp: Optional lead SNP identifier for highlighting.
        gwas_threshold: GWAS significance threshold (default 5e-8).
        eqtl_threshold: eQTL significance threshold (default 1e-5).
        show_correlation: Whether to display Pearson correlation.
        color_by_effect: Whether to color by effect direction agreement.
        gwas_effect_col: Column name for GWAS effect sizes.
        eqtl_effect_col: Column name for eQTL effect sizes.
        h4_posterior: Optional COLOC H4 posterior probability to display.
        figsize: Figure size as (width, height).
    """

    model_config = ConfigDict(frozen=True)

    gwas_p_col: str = Field(default="p_gwas", description="GWAS p-value column")
    eqtl_p_col: str = Field(default="p_eqtl", description="eQTL p-value column")
    pos_col: str = Field(default="pos", description="Position column")
    rs_col: Optional[str] = Field(default="rs", description="SNP ID column")
    ld_col: Optional[str] = Field(default=None, description="Pre-computed LD column")
    lead_snp: Optional[str] = Field(default=None, description="Lead SNP ID")
    gwas_threshold: float = Field(
        default=5e-8, gt=0, le=1, description="GWAS significance"
    )
    eqtl_threshold: float = Field(
        default=1e-5, gt=0, le=1, description="eQTL significance"
    )
    show_correlation: bool = Field(default=True, description="Show Pearson correlation")
    color_by_effect: bool = Field(
        default=False, description="Color by effect agreement"
    )
    gwas_effect_col: Optional[str] = Field(
        default=None, description="GWAS effect column"
    )
    eqtl_effect_col: Optional[str] = Field(
        default=None, description="eQTL effect column"
    )
    h4_posterior: Optional[float] = Field(
        default=None, ge=0, le=1, description="COLOC H4 PP"
    )
    figsize: Tuple[float, float] = Field(default=(8.0, 8.0), description="Figure size")

    @model_validator(mode="after")
    def validate_effect_coloring(self) -> "ColocConfig":
        """Validate that effect coloring has required columns."""
        if self.color_by_effect:
            if self.gwas_effect_col is None or self.eqtl_effect_col is None:
                raise ValueError(
                    "color_by_effect=True requires gwas_effect_col and eqtl_effect_col"
                )
        return self


__all__ = [
    "RegionConfig",
    "ColumnConfig",
    "DisplayConfig",
    "LDConfig",
    "PlotConfig",
    "StackedPlotConfig",
    "ColocConfig",
]

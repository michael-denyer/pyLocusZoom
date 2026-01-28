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

from typing import List, Optional, Tuple

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

        Note: For StackedPlotConfig, ld_reference_file may be provided without
        lead_pos when lead_positions list is used (broadcast mode). This is
        validated at the StackedPlotConfig level, not here.
        """
        # Validation moved to StackedPlotConfig.validate_broadcast_ld
        # to allow broadcast mode where lead_positions list is used instead
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
    def from_kwargs(
        cls,
        *,
        # Region params (required)
        chrom: int,
        start: int,
        end: int,
        # Column params
        pos_col: str = "ps",
        p_col: str = "p_wald",
        rs_col: str = "rs",
        # Display params
        snp_labels: bool = True,
        label_top_n: int = 5,
        show_recombination: bool = True,
        figsize: Tuple[float, float] = (12.0, 8.0),
        # LD params
        lead_pos: Optional[int] = None,
        ld_reference_file: Optional[str] = None,
        ld_col: Optional[str] = None,
    ) -> "PlotConfig":
        """Create PlotConfig from flat keyword arguments.

        Factory method that accepts parameters matching the plot() method
        signature, enabling backward compatibility with existing code.

        Args:
            chrom: Chromosome number.
            start: Start position (bp).
            end: End position (bp).
            pos_col: Column name for position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP ID.
            snp_labels: Whether to show SNP labels.
            label_top_n: Number of top SNPs to label.
            show_recombination: Whether to show recombination overlay.
            figsize: Figure size (width, height).
            lead_pos: Position of lead SNP.
            ld_reference_file: PLINK binary fileset path.
            ld_col: Pre-computed LD column name.

        Returns:
            PlotConfig with nested config objects.

        Raises:
            ValidationError: If parameters are invalid.
        """
        return cls(
            region=RegionConfig(chrom=chrom, start=start, end=end),
            columns=ColumnConfig(pos_col=pos_col, p_col=p_col, rs_col=rs_col),
            display=DisplayConfig(
                snp_labels=snp_labels,
                label_top_n=label_top_n,
                show_recombination=show_recombination,
                figsize=figsize,
            ),
            ld=LDConfig(
                lead_pos=lead_pos,
                ld_reference_file=ld_reference_file,
                ld_col=ld_col,
            ),
        )


class StackedPlotConfig(BaseModel):
    """Composite configuration for plot_stacked() method.

    Extends PlotConfig pattern with list-based parameters for stacked plots.
    Supports multiple lead positions, panel labels, and LD reference files.

    Attributes:
        region: Genomic region specification (required).
        columns: DataFrame column name mappings.
        display: Display and visual options.
        ld: Linkage disequilibrium configuration (single file for broadcast).
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

    model_config = ConfigDict(frozen=True)

    region: RegionConfig
    columns: ColumnConfig = Field(default_factory=ColumnConfig)
    display: DisplayConfig = Field(default_factory=DisplayConfig)
    ld: LDConfig = Field(default_factory=LDConfig)

    # Stacked-specific list parameters
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
    def validate_broadcast_ld(self) -> "StackedPlotConfig":
        """Validate broadcast LD configuration for stacked plots.

        When ld_reference_file is provided for broadcast, lead_positions must
        be provided to specify the reference SNP for each panel.
        """
        if self.ld.ld_reference_file is not None and self.ld.lead_pos is None:
            # Broadcast mode: ld_reference_file without lead_pos in LDConfig
            # Requires lead_positions list instead
            if self.lead_positions is None:
                raise ValueError(
                    "lead_positions is required when ld_reference_file is provided "
                    "for broadcast (one lead position per panel)"
                )
        return self

    @classmethod
    def from_kwargs(
        cls,
        *,
        # Region params (required)
        chrom: int,
        start: int,
        end: int,
        # Column params
        pos_col: str = "ps",
        p_col: str = "p_wald",
        rs_col: str = "rs",
        # Display params
        snp_labels: bool = True,
        label_top_n: int = 3,  # Default for stacked is 3 (less crowded)
        show_recombination: bool = True,
        figsize: Tuple[float, float] = (12.0, 8.0),
        # LD params (single for broadcast)
        ld_reference_file: Optional[str] = None,
        ld_col: Optional[str] = None,
        # Stacked-specific list params
        lead_positions: Optional[List[int]] = None,
        panel_labels: Optional[List[str]] = None,
        ld_reference_files: Optional[List[str]] = None,
    ) -> "StackedPlotConfig":
        """Create StackedPlotConfig from flat keyword arguments.

        Factory method that accepts parameters matching the plot_stacked()
        method signature, enabling backward compatibility.

        Args:
            chrom: Chromosome number.
            start: Start position (bp).
            end: End position (bp).
            pos_col: Column name for position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP ID.
            snp_labels: Whether to show SNP labels.
            label_top_n: Number of top SNPs to label (default 3 for stacked).
            show_recombination: Whether to show recombination overlay.
            figsize: Figure size (width, height).
            ld_reference_file: Single PLINK fileset (broadcast to all panels).
            ld_col: Pre-computed LD column name.
            lead_positions: List of lead SNP positions.
            panel_labels: List of panel labels.
            ld_reference_files: List of PLINK filesets.

        Returns:
            StackedPlotConfig with nested config objects.

        Raises:
            ValidationError: If parameters are invalid.
        """
        return cls(
            region=RegionConfig(chrom=chrom, start=start, end=end),
            columns=ColumnConfig(pos_col=pos_col, p_col=p_col, rs_col=rs_col),
            display=DisplayConfig(
                snp_labels=snp_labels,
                label_top_n=label_top_n,
                show_recombination=show_recombination,
                figsize=figsize,
            ),
            ld=LDConfig(
                ld_reference_file=ld_reference_file,
                ld_col=ld_col,
            ),
            lead_positions=lead_positions,
            panel_labels=panel_labels,
            ld_reference_files=ld_reference_files,
        )


__all__ = [
    "RegionConfig",
    "ColumnConfig",
    "DisplayConfig",
    "LDConfig",
    "PlotConfig",
    "StackedPlotConfig",
]

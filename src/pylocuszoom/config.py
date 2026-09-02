"""Pydantic configuration classes for pyLocusZoom plot methods.

Each option a plot method accepts is declared once, on the model that owns
it, and the method takes the model as a value. Every model is immutable
(frozen), so one built in a notebook can be handed to as many calls as
needed.

Example:
    >>> from pylocuszoom import DisplayConfig, LDConfig, LocusZoomPlotter
    >>> plotter = LocusZoomPlotter(species="canine")
    >>> fig = plotter.plot(
    ...     gwas_df,
    ...     chrom=1,
    ...     start=1000000,
    ...     end=2000000,
    ...     ld=LDConfig(lead_pos=1500000, ld_col="R2"),
    ...     display=DisplayConfig(snp_labels=False),
    ... )
"""

from typing import Annotated, ClassVar, List, Optional, Tuple, Union

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from ._plotter_utils import DEFAULT_EQTL_THRESHOLD, DEFAULT_GENOMEWIDE_THRESHOLD

PValueThreshold = Annotated[float, Field(gt=0, le=1)]


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
        label_top_n: Number of top SNPs to label per association panel.
            None takes the method's default: 5 on ``plot()``, 3 on
            ``plot_stacked()``, whose panels are shorter.
        show_recombination: Whether to show recombination rate overlay.
        auto_genes: Fetch the gene track when no ``genes_df`` is supplied.
            None inherits the plotter's constructor setting.
        figsize: Figure size as (width, height) in inches.
    """

    model_config = ConfigDict(frozen=True)

    snp_labels: bool = Field(default=True, description="Show SNP labels")
    label_top_n: Optional[int] = Field(
        default=None, ge=0, description="Number of top SNPs to label"
    )
    show_recombination: bool = Field(
        default=True, description="Show recombination overlay"
    )
    auto_genes: Optional[bool] = Field(
        default=None, description="Fetch the gene track; None inherits the plotter"
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


class PanelInputs(BaseModel):
    """Caller-supplied data for the optional panels beneath the association track."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    genes_df: Optional[pd.DataFrame] = Field(
        default=None, description="Gene annotations"
    )
    exons_df: Optional[pd.DataFrame] = Field(
        default=None, description="Exon annotations"
    )
    recomb_df: Optional[pd.DataFrame] = Field(
        default=None, description="Recombination rates"
    )
    eqtl_df: Optional[pd.DataFrame] = Field(default=None, description="eQTL results")
    eqtl_gene: Optional[str] = Field(default=None, description="eQTL gene filter")
    eqtl_threshold: float = Field(
        default=DEFAULT_EQTL_THRESHOLD, description="eQTL significance"
    )
    finemapping_df: Optional[pd.DataFrame] = Field(
        default=None, description="Fine-mapping results"
    )
    finemapping_cs_col: Optional[str] = Field(
        default="cs", description="Credible-set column"
    )
    ld_heatmap_df: Optional[pd.DataFrame] = Field(default=None, description="LD matrix")
    ld_heatmap_snp_ids: Optional[List[str]] = Field(
        default=None, description="LD matrix row/column SNP ids"
    )
    ld_heatmap_height: float = Field(
        default=0.25,
        description="Heatmap height as a fraction of the association panel",
    )
    ld_heatmap_metric: str = Field(default="r2", description="LD metric label")

    @model_validator(mode="after")
    def validate_heatmap_requires_snp_ids(self) -> "PanelInputs":
        """Validate that an LD heatmap matrix names its rows and columns."""
        if self.ld_heatmap_df is not None and self.ld_heatmap_snp_ids is None:
            raise ValueError(
                "ld_heatmap_snp_ids is required when ld_heatmap_df is provided"
            )
        return self


class PlotConfig(BaseModel):
    """Everything ``plot()`` was asked for, as one validated value.

    ``plot()`` builds one from its arguments; the cross-model rules live
    here. Callers do not construct it.

    Attributes:
        region: Genomic region specification (required).
        columns: DataFrame column name mappings.
        display: Display and visual options.
        ld: Linkage disequilibrium configuration.
        panels: Data for the optional panels beneath the association track.
    """

    model_config = ConfigDict(frozen=True)

    region: RegionConfig
    columns: ColumnConfig = Field(default_factory=ColumnConfig)
    display: DisplayConfig = Field(default_factory=DisplayConfig)
    ld: LDConfig = Field(default_factory=LDConfig)
    panels: PanelInputs = Field(default_factory=PanelInputs)

    _lead_required_error: ClassVar[str] = (
        "lead_pos is required when ld_reference_file is provided"
    )

    def _lead_is_set(self) -> bool:
        """Report whether a lead position is available to compute LD against."""
        return self.ld.lead_pos is not None

    @model_validator(mode="after")
    def validate_ld_requires_lead(self) -> "PlotConfig":
        """Validate that computing LD from a fileset has a lead to compute against."""
        if self.ld.ld_reference_file is not None and not self._lead_is_set():
            raise ValueError(self._lead_required_error)
        return self


class StackedPlotConfig(PlotConfig):
    """Everything ``plot_stacked()`` was asked for, as one validated value.

    Extends :class:`PlotConfig` with the per-panel lists, each of which must
    hold one entry per GWAS frame.

    Attributes:
        n_panels: Number of association panels, one per GWAS frame.
        lead_positions: List of lead SNP positions (one per panel).
        panel_labels: List of panel labels (one per panel).
        ld_reference_files: List of PLINK filesets (one per panel).
    """

    n_panels: int = Field(..., ge=1, description="Number of association panels")
    lead_positions: Optional[List[int]] = Field(
        default=None, description="Lead SNP positions (one per panel)"
    )
    panel_labels: Optional[List[str]] = Field(
        default=None, description="Panel labels (one per panel)"
    )
    ld_reference_files: Optional[List[str]] = Field(
        default=None, description="PLINK filesets (one per panel)"
    )

    _lead_required_error: ClassVar[str] = (
        "lead_positions is required when ld_reference_file is provided "
        "for broadcast (one lead position per panel)"
    )

    def _lead_is_set(self) -> bool:
        """A stacked figure may take one lead per panel instead of one overall."""
        return self.ld.lead_pos is not None or self.lead_positions is not None

    @model_validator(mode="after")
    def validate_per_panel_lists(self) -> "StackedPlotConfig":
        """Validate that each per-panel list has one entry per panel."""
        for name in ("lead_positions", "panel_labels", "ld_reference_files"):
            value = getattr(self, name)
            if value is not None and len(value) != self.n_panels:
                raise ValueError(
                    f"{name} length ({len(value)}) must match "
                    f"number of GWAS DataFrames ({self.n_panels})"
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
        gwas_threshold: GWAS significance threshold, or None to draw no line.
        eqtl_threshold: eQTL significance threshold, or None to draw no line.
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
    gwas_threshold: Optional[PValueThreshold] = Field(
        default=DEFAULT_GENOMEWIDE_THRESHOLD, description="GWAS significance"
    )
    eqtl_threshold: Optional[PValueThreshold] = Field(
        default=DEFAULT_EQTL_THRESHOLD, description="eQTL significance"
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
    "PanelInputs",
    "PlotConfig",
    "StackedPlotConfig",
    "ColocConfig",
]

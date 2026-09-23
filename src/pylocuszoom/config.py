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

import os
import warnings
from typing import (
    Annotated,
    Any,
    ClassVar,
    List,
    Literal,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

import matplotlib.colors as mcolors
import pandas as pd
import pydantic
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from ._liftover import CoordinateLifter, load_chain
from ._plotter_utils import (
    CHROMOSOME_GAP,
    DEFAULT_EQTL_THRESHOLD,
    DEFAULT_GENOMEWIDE_THRESHOLD,
)
from .exceptions import ValidationError
from .schemas import (
    DEPRECATED_ALIAS_REMOVED_IN,
    DEPRECATED_COLUMN_ALIASES,
    Canonical,
)

PValueThreshold = Annotated[float, Field(gt=0, le=1)]


def _describe(error: pydantic.ValidationError) -> str:
    """Name each failing field and why, without pydantic's type codes and URLs."""
    faults = []
    for fault in error.errors():
        field = ".".join(str(part) for part in fault["loc"])
        message = fault["msg"].removeprefix("Value error, ")
        faults.append(f"{field}: {message}" if field else message)
    return f"Invalid {error.title}: " + "; ".join(faults)


class _Config(BaseModel):
    """Base of every config model: a rejected value raises our ValidationError.

    pydantic's own ``ValidationError`` shares the name but not the hierarchy,
    so ``except PyLocusZoomError`` would miss it.
    """

    def __init__(self, /, **data: Any) -> None:
        try:
            super().__init__(**data)
        except pydantic.ValidationError as error:
            raise ValidationError(_describe(error)) from error

    def __setattr__(self, name: str, value: Any) -> None:
        try:
            super().__setattr__(name, value)
        except pydantic.ValidationError as error:
            raise ValidationError(_describe(error)) from error


class RegionConfig(_Config):
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


class ColumnConfig(_Config):
    """DataFrame column name mappings for GWAS data.

    The defaults are the canonical names every loader emits, so a loaded
    frame needs no config at all.

    Attributes:
        pos_col: Column name for genomic position.
        p_col: Column name for p-value.
        rs_col: Column name for SNP identifier.
    """

    model_config = ConfigDict(frozen=True)

    pos_col: str = Field(default=Canonical.POS, description="Position column name")
    p_col: str = Field(default=Canonical.P, description="P-value column name")
    rs_col: str = Field(default=Canonical.RS, description="SNP ID column name")


class DisplayConfig(_Config):
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

    def with_defaults(self, *, label_top_n: int, auto_genes: bool) -> "DisplayConfig":
        """Return a copy with the unset fields filled from the caller's defaults.

        Args:
            label_top_n: The plot method's own label count.
            auto_genes: The plotter's constructor setting.
        """
        return self.model_copy(
            update={
                "label_top_n": (
                    label_top_n if self.label_top_n is None else self.label_top_n
                ),
                "auto_genes": auto_genes
                if self.auto_genes is None
                else self.auto_genes,
            }
        )


class LDConfig(_Config):
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


class LiftoverConfig(_Config):
    """Plot summary statistics from one genome build on another build's annotations.

    With a ``lifter`` or a ``chain_path``, ``plot()`` treats ``gwas_df``,
    ``start``, ``end`` and ``ld.lead_pos`` as source-build coordinates and
    lifts the region's SNPs to the plotter's ``genome_build`` before drawing.
    Gene, exon, eQTL and fine-mapping frames are taken as already in the
    target build. The default, neither set, plots without liftover.

    Attributes:
        lifter: Coordinate lifter exposing pyliftover's ``convert_coordinate``,
            such as a ``pyliftover.LiftOver``.
        chain_path: UCSC chain file from the source build to the plotter's
            build, loaded with pyliftover. Loaded chains are cached per path.
        lift_recombination: Lift the recombination maps the plotter loads,
            managed or from ``recomb_data_dir``, through the same chain rather
            than the registered one. The maps must be in the chain's source
            build. A ``PanelInputs.recomb_df`` is never lifted.
    """

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    lifter: Optional[CoordinateLifter] = Field(
        default=None, description="Coordinate lifter"
    )
    chain_path: Optional[Union[str, os.PathLike]] = Field(
        default=None, description="UCSC chain file"
    )
    lift_recombination: bool = Field(
        default=False, description="Lift the recombination maps too"
    )

    @model_validator(mode="after")
    def validate_one_source(self) -> "LiftoverConfig":
        """Validate that at most one lifter source is named, and one is if needed."""
        if self.lifter is not None and self.chain_path is not None:
            raise ValueError("Pass either lifter or chain_path, not both")
        if self.lift_recombination and self.lifter is None and self.chain_path is None:
            raise ValueError("lift_recombination requires a lifter or chain_path")
        return self

    def resolve(self) -> Optional[CoordinateLifter]:
        """Return the lifter to use, loading ``chain_path`` if needed, or None."""
        if self.chain_path is not None:
            return load_chain(self.chain_path)
        return self.lifter


class PanelInputs(_Config):
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


class PlotConfig(_Config):
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
    lead_positions: Optional[List[Annotated[int, Field(ge=1)]]] = Field(
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


class GenomeWideConfig(_Config):
    """Column names and chromosome order for the genome-wide plot families.

    Manhattan, QQ and Miami plots lay a whole-genome frame out along one
    chromosome axis, so they share one column contract. ``plot_qq`` reads
    only ``p_col``.

    Attributes:
        chrom_col: Column name for chromosome.
        pos_col: Column name for genomic position.
        p_col: Column name for p-value.
        custom_chrom_order: Chromosome order along the axis, overriding the
            plotter's species order.
    """

    model_config = ConfigDict(frozen=True)

    chrom_col: str = Field(
        default=Canonical.CHROM, description="Chromosome column name"
    )
    pos_col: str = Field(default=Canonical.POS, description="Position column name")
    p_col: str = Field(default=Canonical.P, description="P-value column name")
    custom_chrom_order: Optional[List[str]] = Field(
        default=None, description="Chromosome order overriding the species"
    )


PositiveFontSize = Annotated[int, Field(gt=0)]


class GenomeWideStyle(_Config):
    """Colours, points, fonts and chromosome axis of the genome-wide plots.

    Every Manhattan, QQ, Manhattan-QQ, stacked and Miami method takes one as
    ``style``. The default reproduces each method's own look, so a field left
    unset changes nothing. A field whose default is None takes the method's
    own value, which differs by method: the chromosome ticks
    are 8 pt on a genomic axis and 10 pt on a category axis, the points of a
    categorical Manhattan are larger, and a stacked figure uses smaller panel
    titles and axis labels.

    Attributes:
        palette: Colours cycled over the chromosomes in display order, or
            over the categories of a categorical Manhattan. Any matplotlib
            colour spec; each is stored as a hex string, which every backend
            accepts. To use a matplotlib colormap, pass its ``colors``. None
            keeps the default glasbey palette.
        point_size: Marker area in matplotlib's ``s`` units, applied to the
            Manhattan and QQ points. The interactive backends convert it to
            a diameter. None keeps the method's size.
        point_alpha: Marker opacity in (0, 1] for Manhattan and QQ points.
            None draws them opaque.
        title_fontsize: Size of the figure title that ``title`` sets on a
            multi-panel figure (Manhattan-QQ, stacked, Miami).
        panel_title_fontsize: Size of the title drawn on each panel: the
            "Manhattan Plot" and QQ lambda titles, and the ``title`` of a
            single-panel ``plot_manhattan`` or ``plot_qq``.
        axis_label_fontsize: Size of the x and y axis labels.
        tick_label_fontsize: Size of the tick labels on both axes.
        tick_step: Label every ``tick_step``-th chromosome (or category) that
            carries data, starting with the first. At least 1.
        tick_rotation: Rotation of the chromosome or category tick labels in
            degrees. None keeps the method's rotation.
        chrom_gap: Gap in base pairs between one chromosome's last position
            and the next chromosome's first on a genomic axis.
        line_style: Matplotlib linestyle of the significance and suggestive
            lines and the QQ diagonal: ``"-"``, ``"--"``, ``":"`` or ``"-."``.
        line_width: Width of the same lines.
        title_fontweight: Weight of the figure title and the panel titles,
            ``"bold"`` or ``"normal"``.
        point_edge_width: Outline width of the Manhattan and QQ points; 0
            draws no outline. None keeps the method's width.
        y_headroom: Space left above the highest Manhattan point or threshold
            line, as a fraction of its height.
        manhattan_qq_width_ratio: Width of the Manhattan panel relative to the
            QQ panel in a Manhattan-QQ figure.
    """

    model_config = ConfigDict(frozen=True)

    palette: Optional[Tuple[str, ...]] = Field(
        default=None, description="Colours cycled over chromosomes"
    )
    point_size: Optional[float] = Field(default=None, gt=0, description="Marker area")
    point_alpha: Optional[float] = Field(
        default=None, gt=0, le=1, description="Marker opacity"
    )
    title_fontsize: PositiveFontSize = Field(
        default=14, description="Figure title size"
    )
    panel_title_fontsize: Optional[PositiveFontSize] = Field(
        default=None, description="Panel title size"
    )
    axis_label_fontsize: Optional[PositiveFontSize] = Field(
        default=None, description="Axis label size"
    )
    tick_label_fontsize: Optional[PositiveFontSize] = Field(
        default=None, description="Tick label size"
    )
    tick_step: int = Field(default=1, ge=1, description="Label every n-th tick")
    tick_rotation: Optional[int] = Field(
        default=None, description="Chromosome tick label rotation"
    )
    chrom_gap: int = Field(
        default=CHROMOSOME_GAP, ge=0, description="Gap between chromosomes (bp)"
    )
    line_style: Literal["-", "--", ":", "-."] = Field(
        default="--", description="Threshold line and QQ diagonal style"
    )
    line_width: float = Field(
        default=1.0, gt=0, description="Threshold line and QQ diagonal width"
    )
    title_fontweight: Literal["bold", "normal"] = Field(
        default="bold", description="Figure and panel title weight"
    )
    point_edge_width: Optional[float] = Field(
        default=None, ge=0, description="Marker outline width"
    )
    y_headroom: float = Field(
        default=0.1, ge=0, description="Manhattan space above the top point or line"
    )
    manhattan_qq_width_ratio: float = Field(
        default=2.5, gt=0, description="Manhattan panel width over QQ panel width"
    )

    @field_validator("palette", mode="before")
    @classmethod
    def validate_palette(cls, v: Any) -> Any:
        """Validate a non-empty sequence of colours and store each as hex."""
        if v is None:
            return v
        if isinstance(v, str):
            raise ValueError("palette must be a sequence of colours, not one string")
        colours = list(v)
        if not colours:
            raise ValueError("palette must name at least one colour")
        for colour in colours:
            if not mcolors.is_color_like(colour):
                raise ValueError(f"palette entry {colour!r} is not a colour")
        return tuple(mcolors.to_hex(colour) for colour in colours)


class ColocConfig(_Config):
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
    pos_col: str = Field(default=Canonical.POS, description="Position column")
    rs_col: Optional[str] = Field(default=Canonical.RS, description="SNP ID column")
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


Columns = TypeVar("Columns", ColumnConfig, GenomeWideConfig)


def resolve_deprecated_columns(
    df: pd.DataFrame,
    columns: Columns,
    *,
    fields: tuple[str, ...] = ("pos_col", "p_col"),
) -> Columns:
    """Fall back to a frame's pre-4.0 column names, with a deprecation warning.

    A GWAS frame written by a 3.x loader carries ``ps`` and ``p_wald`` where a
    4.0 one carries the canonical ``pos`` and ``p_value``. Where the canonical
    column is absent and the old spelling is present, the old spelling is used
    and a ``DeprecationWarning`` names its replacement. A caller who asked for
    a column name of their own is left alone, because only the canonical names
    have aliases.

    Args:
        df: The frame about to be plotted.
        columns: The column model the plot method was given.
        fields: Column roles consumed by this input boundary.

    Returns:
        ``columns`` unchanged, or a copy naming the frame's old columns.
    """
    updates = {}
    for field in fields:
        name = getattr(columns, field)
        alias = DEPRECATED_COLUMN_ALIASES.get(name)
        if alias is None or name in df.columns or alias not in df.columns:
            continue
        warnings.warn(
            f"Column '{alias}' is the pre-4.0 name for '{name}'. pyLocusZoom "
            f"loaders now emit '{name}'; rename the column, or pass the name "
            f"you want. The fallback is removed in "
            f"{DEPRECATED_ALIAS_REMOVED_IN}.",
            DeprecationWarning,
            stacklevel=3,
        )
        updates[field] = alias
    return columns.model_copy(update=updates) if updates else columns


__all__ = [
    "RegionConfig",
    "ColumnConfig",
    "DisplayConfig",
    "LDConfig",
    "LiftoverConfig",
    "PanelInputs",
    "PlotConfig",
    "StackedPlotConfig",
    "GenomeWideConfig",
    "GenomeWideStyle",
    "ColocConfig",
    "resolve_deprecated_columns",
]

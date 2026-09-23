"""Main LocusZoomPlotter class for regional association plots.

Orchestrates all components (LD coloring, gene track, recombination overlay,
SNP labels) into a unified plotting interface.

Supports multiple backends:
- matplotlib (default): Static publication-quality plots
- plotly: Interactive HTML with hover tooltips
- bokeh: Interactive HTML for dashboards
"""

import warnings
from dataclasses import dataclass
from typing import Any, List, Optional, Union

import pandas as pd

from ._data import prepare_pvalue_data
from ._figure import FigurePlan, render_figure
from ._ld_plotting import enrich_with_ld
from ._liftover import CoordinateLifter, lift_window
from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from .backends import BackendType, get_backend
from .config import (
    ColumnConfig,
    DisplayConfig,
    LDConfig,
    LiftoverConfig,
    PanelInputs,
    PlotConfig,
    RegionConfig,
    StackedPlotConfig,
)
from .exceptions import ReferenceAPIError, ValidationError
from .ld import find_plink
from .logging import enable_logging, logger
from .panels import (
    AssociationPanel,
    EqtlPanel,
    FinemappingPanel,
    GenePanel,
    HeatmapPanel,
    RegionalPanel,
    hover_for_association,
)
from .recombination import RecombResult, RecombStatus, recomb_for_region
from .reference_genes import get_genes_for_build, source_for
from .schemas import Canonical, gwas_plot_spec
from .species import Species, resolve_species
from .utils import DataFrameLike, filter_by_region, to_pandas
from .validation import check, resolve_column


@dataclass(frozen=True)
class _AssociationInput:
    """One region-selected frame with its effective per-panel options."""

    data: pd.DataFrame
    columns: ColumnConfig
    rs_col: Optional[str]
    ld: LDConfig
    lead_index: Optional[int]
    label: Optional[str] = None

    @classmethod
    def prepare(
        cls,
        frame: pd.DataFrame,
        region: RegionConfig,
        columns: ColumnConfig,
        ld: LDConfig,
        label: Optional[str] = None,
    ) -> "_AssociationInput":
        check(frame, gwas_plot_spec(columns.pos_col, columns.p_col))
        resolve_column(frame, ld.ld_col, parameter="ld_col")
        rs_col = resolve_column(
            frame, columns.rs_col, parameter="rs_col", optional_default=Canonical.RS
        )
        if ld.ld_reference_file is not None and rs_col is None:
            raise ValidationError(
                "ld_reference_file needs SNP ids to compute LD, and column "
                f"'{columns.rs_col}' is not in the GWAS data. Add it, or name "
                "the SNP id column with ColumnConfig(rs_col=...)."
            )
        selected = filter_by_region(
            frame,
            region=(region.chrom, region.start, region.end),
            chrom_col=columns.chrom_col,
            pos_col=columns.pos_col,
        )
        data = prepare_pvalue_data(selected, columns.p_col, "regional")
        data = data.reset_index(drop=True)
        candidates = (
            data if ld.lead_pos is None else data[data[columns.pos_col] == ld.lead_pos]
        )
        lead_index = (
            int(candidates["neglog10p"].idxmax()) if not candidates.empty else None
        )
        if ld.lead_pos is not None and lead_index is None:
            logger.warning(
                "Lead SNP at position {} not found in region; LD coloring will be skipped",
                ld.lead_pos,
            )
        return cls(data, columns, rs_col, ld, lead_index, label)


class LocusZoomPlotter:
    """Regional association plot generator with LD coloring and annotations.

    Creates LocusZoom-style regional plots with:
    - LD coloring based on R-squared with lead variant
    - Gene and exon tracks
    - Recombination rate overlays (canine built-in, or user-provided)
    - Automatic SNP labeling

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        species: Species name, alias or record ('canine', 'dog', 'feline',
            'human', ..., or None for custom). An unknown name raises
            ValidationError. Canine has built-in recombination maps.
        genome_build: Genome build for coordinate system. For canine:
            "canfam3.1" (default) or "canfam4". If "canfam4", recombination
            maps are automatically lifted over from CanFam3.1.
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').
            Defaults to 'matplotlib' for static plots.
        plink_path: Path to PLINK executable for LD calculation.
            Auto-detects if None.
        recomb_data_dir: Directory containing recombination maps.
            Uses platform cache if None.
        genomewide_threshold: P-value threshold for significance line.
        log_level: Logging level ("DEBUG", "INFO", "WARNING", "ERROR", or None
            to disable). Defaults to "INFO".

    Example:
        >>> # Static plot (default)
        >>> plotter = LocusZoomPlotter(species="canine")
        >>>
        >>> # Interactive plot with plotly
        >>> plotter = LocusZoomPlotter(species="canine", backend="plotly")
        >>>
        >>> from pylocuszoom import LDConfig
        >>> fig = plotter.plot(
        ...     gwas_df,
        ...     chrom=1,
        ...     start=1000000,
        ...     end=2000000,
        ...     ld=LDConfig(lead_pos=1500000),
        ... )
        >>> fig.savefig("regional_plot.png", dpi=150)  # matplotlib
        >>> # or fig.save("plot.html")  # plotly/bokeh
    """

    def __init__(
        self,
        species: str | Species | None = "canine",
        genome_build: Optional[str] = None,
        backend: BackendType = "matplotlib",
        plink_path: Optional[str] = None,
        recomb_data_dir: Optional[str] = None,
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
        log_level: Optional[str] = "INFO",
        auto_genes: bool = False,
    ):
        """Initialize the plotter."""
        if log_level is not None:
            enable_logging(log_level)

        self.species = resolve_species(species)
        self.genome_build = genome_build or (
            self.species.default_build if self.species else None
        )
        self._backend = get_backend(backend)
        self.plink_path = plink_path or find_plink()
        self.recomb_data_dir = recomb_data_dir
        self.genomewide_threshold = genomewide_threshold
        self._auto_genes = auto_genes
        self._recomb_cache = {}

    def _get_recomb_for_region(
        self,
        chrom: int,
        start: int,
        end: int,
        lifter: Optional[CoordinateLifter] = None,
    ) -> RecombResult:
        """Get a region's recombination rates, or the reason there are none.

        Caches per region, build and lifter. The caller renders the outcome;
        this does not warn, so a region asked for twice is reported once.
        """
        cache_key = (chrom, start, end, self.genome_build, lifter)
        if cache_key not in self._recomb_cache:
            self._recomb_cache[cache_key] = recomb_for_region(
                chrom=chrom,
                start=start,
                end=end,
                species=self.species,
                data_dir=self.recomb_data_dir,
                genome_build=self.genome_build,
                lifter=lifter,
            )
        return self._recomb_cache[cache_key]

    def plot(
        self,
        gwas_df: DataFrameLike,
        *,
        chrom: Union[int, str],
        start: int,
        end: int,
        columns: ColumnConfig = ColumnConfig(),
        display: DisplayConfig = DisplayConfig(),
        ld: LDConfig = LDConfig(),
        panels: PanelInputs = PanelInputs(),
        significance_threshold: ThresholdArg = UNSET,
        liftover: LiftoverConfig = LiftoverConfig(),
    ) -> Any:
        """Create a regional association plot for a single locus.

        Plots ``-log10(p)`` against genomic position for the specified region,
        optionally overlaid with LD colouring, recombination rate, SNP labels,
        a gene track, and fine-mapping, eQTL, or LD-heatmap panels beneath.
        Every option is declared once, on the config model that owns it; the
        models are frozen, so one built in a notebook serves many calls.

        Args:
            gwas_df: GWAS summary statistics. Must contain the columns named
                by ``columns.pos_col`` and ``columns.p_col``; needs
                ``columns.rs_col`` only when LD is computed from a fileset.
            chrom: Chromosome of the region.
            start: Region start position (bp, inclusive, ``>= 1``).
            end: Region end position (bp, inclusive, ``> start``).
            columns: :class:`~pylocuszoom.ColumnConfig` naming the position,
                p-value and SNP id columns of ``gwas_df``.
            display: :class:`~pylocuszoom.DisplayConfig` for SNP labels, the
                recombination overlay, automatic gene fetching and the figure
                size. ``label_top_n`` defaults to 5 here.
            ld: :class:`~pylocuszoom.LDConfig` naming the lead SNP and the LD
                source. ``lead_pos`` is auto-detected as the strongest
                in-region p-value when omitted, and is required when
                ``ld_reference_file`` is set and ``ld_col`` is not.
            panels: :class:`~pylocuszoom.PanelInputs` carrying the frames for
                the optional gene, eQTL, fine-mapping and LD-heatmap panels,
                plus a caller-supplied recombination frame.
            significance_threshold: P-value for the genome-wide significance
                line. Defaults to the plotter's ``genomewide_threshold``;
                pass None to draw no line.
            liftover: :class:`~pylocuszoom.LiftoverConfig` naming a chain when
                ``gwas_df``, ``start``, ``end`` and ``ld.lead_pos`` are in
                another build than the plotter's ``genome_build``. The region's
                SNPs are lifted with :func:`~pylocuszoom.liftover_region` and
                the window keeps the requested margins around the outermost
                lifted SNPs.

        Returns:
            Backend-specific figure object (``matplotlib.figure.Figure``,
            ``plotly.graph_objects.Figure``, or ``bokeh.layouts.Column``).

        Raises:
            ValidationError: On an invalid region or a contradictory config,
                a missing required GWAS column, or when no SNP in the region
                lifts to the target build.
            pylocuszoom.exceptions.PlinkError: When PLINK itself fails
                (timeout, non-zero exit, corrupt ``.bed``, missing output).
                The specific "empty LD output" case, a singleton lead SNP with
                no neighbours in the window, is downgraded to a warning and
                the plot is drawn without LD colouring.

        Example:
            >>> from pylocuszoom import LDConfig, PanelInputs
            >>> fig = plotter.plot(
            ...     gwas_df,
            ...     chrom=1,
            ...     start=1_000_000,
            ...     end=2_000_000,
            ...     ld=LDConfig(lead_pos=1_500_000, ld_reference_file="ref"),
            ...     panels=PanelInputs(genes_df=genes_df),
            ... )
        """
        gwas_df = to_pandas(gwas_df)
        config = PlotConfig(
            region=RegionConfig(chrom=chrom, start=start, end=end),
            columns=columns,
            display=display,
            ld=ld,
            panels=panels,
        )
        lifter = liftover.resolve()
        if lifter is not None:
            [gwas_df], region, [ld] = self._lift(
                [gwas_df], config.region, columns, [ld], lifter
            )
            config = config.model_copy(update={"region": region, "ld": ld})
        return self._render_regional(
            config,
            [_AssociationInput.prepare(gwas_df, config.region, columns, ld)],
            threshold=resolve_threshold(
                significance_threshold, self.genomewide_threshold
            ),
            label_top_n=5,
            association_height=display.figsize[1] * 0.6,
            min_figure_height=0.0,
            recomb_lifter=lifter if liftover.lift_recombination else None,
        )

    def _lift(
        self,
        frames: List[pd.DataFrame],
        region: RegionConfig,
        columns: ColumnConfig,
        lds: List[LDConfig],
        lifter: CoordinateLifter,
    ) -> tuple[List[pd.DataFrame], RegionConfig, List[LDConfig]]:
        """Lift validated source-build panels to the plotter's build.

        Returns the lifted frames, the lifted window and each panel's LD
        config with its lead lifted. A lead that does not lift is
        auto-detected instead, unless the panel computes LD from a fileset,
        which needs the lead; that is an error rather than a warning the
        next check would contradict.
        """
        window = lift_window(
            frames,
            chrom=region.chrom,
            start=region.start,
            end=region.end,
            chrom_col=columns.chrom_col,
            pos_col=columns.pos_col,
            lead_positions=[ld.lead_pos for ld in lds],
            lifter=lifter,
            build=self.genome_build,
        )
        for ld, lead in zip(lds, window.lead_positions):
            if ld.lead_pos is not None and lead is None and ld.ld_reference_file:
                raise ValidationError(
                    f"Lead SNP at chr{region.chrom}:{ld.lead_pos} did not lift to "
                    f"{self.genome_build}, and LD from ld_reference_file needs a "
                    "lead; pass a lead that lifts or drop ld_reference_file"
                )
        for note in window.notes:
            warnings.warn(note, stacklevel=3)
        return (
            window.frames,
            RegionConfig(chrom=region.chrom, start=window.start, end=window.end),
            [
                ld.model_copy(update={"lead_pos": lead})
                for ld, lead in zip(lds, window.lead_positions)
            ],
        )

    def plot_stacked(
        self,
        gwas_dfs: List[DataFrameLike],
        *,
        chrom: Union[int, str],
        start: int,
        end: int,
        columns: ColumnConfig = ColumnConfig(),
        display: DisplayConfig = DisplayConfig(),
        ld: LDConfig = LDConfig(),
        panels: PanelInputs = PanelInputs(),
        lead_positions: Optional[List[int]] = None,
        panel_labels: Optional[List[str]] = None,
        ld_reference_files: Optional[List[str]] = None,
        significance_threshold: ThresholdArg = UNSET,
        liftover: LiftoverConfig = LiftoverConfig(),
    ) -> Any:
        """Create stacked regional association plots for multiple GWAS.

        Each frame in ``gwas_dfs`` becomes one association panel; optional
        fine-mapping, eQTL, gene-track, and LD-heatmap panels follow beneath.
        The config models behave as documented on :meth:`plot`, with two
        per-panel differences: ``display.label_top_n`` defaults to 3, and
        ``ld.lead_pos`` and ``ld.ld_reference_file`` apply to every panel
        unless the per-panel lists below override them.

        Args:
            gwas_dfs: One GWAS summary-statistics frame per panel.
            lead_positions: One lead position per panel. Auto-detected as
                the strongest in-region p-value when omitted. Every panel that
                computes LD from a reference fileset needs a lead, from this
                list or from ``ld.lead_pos``.
            panel_labels: One label per panel, or None for none.
            ld_reference_files: One PLINK fileset per panel, replacing the
                broadcast ``ld.ld_reference_file``.
            significance_threshold: As on :meth:`plot`.
            liftover: As on :meth:`plot`, applied to every frame. The window
                spans the lifted SNPs of all panels, with the requested
                margins around them.

        Raises:
            ValidationError: If ``gwas_dfs`` is empty, a per-panel list has a
                different length, or a panel computes LD without a lead.

        Example:
            >>> fig = plotter.plot_stacked(
            ...     [gwas_a, gwas_b],
            ...     chrom=1,
            ...     start=1_000_000,
            ...     end=2_000_000,
            ...     lead_positions=[1_500_000, 1_700_000],
            ...     panel_labels=["Height", "BMI"],
            ...     panels=PanelInputs(genes_df=genes_df),
            ... )
        """
        gwas_dfs = [to_pandas(df) for df in gwas_dfs]
        if not gwas_dfs:
            raise ValidationError("At least one GWAS DataFrame required")
        config = StackedPlotConfig(
            region=RegionConfig(chrom=chrom, start=start, end=end),
            columns=columns,
            display=display,
            ld=ld,
            panels=panels,
            n_panels=len(gwas_dfs),
            lead_positions=lead_positions,
            panel_labels=panel_labels,
            ld_reference_files=ld_reference_files,
        )
        panel_lds = config.panel_lds()
        lifter = liftover.resolve()
        if lifter is not None:
            gwas_dfs, region, panel_lds = self._lift(
                gwas_dfs, config.region, columns, panel_lds, lifter
            )
            config = config.model_copy(update={"region": region})
        association = [
            _AssociationInput.prepare(
                frame,
                config.region,
                columns,
                panel_ld,
                panel_labels[index] if panel_labels is not None else None,
            )
            for index, (frame, panel_ld) in enumerate(zip(gwas_dfs, panel_lds))
        ]
        return self._render_regional(
            config,
            association,
            threshold=resolve_threshold(
                significance_threshold, self.genomewide_threshold
            ),
            label_top_n=3,
            association_height=2.5,
            min_figure_height=display.figsize[1],
            recomb_lifter=lifter if liftover.lift_recombination else None,
        )

    def _render_regional(
        self,
        config: PlotConfig,
        association_inputs: List[_AssociationInput],
        *,
        threshold: Optional[float],
        label_top_n: int,
        association_height: float,
        min_figure_height: float,
        recomb_lifter: Optional[CoordinateLifter] = None,
    ) -> Any:
        """Build the panel plan for one regional figure and render it.

        ``plot()`` and ``plot_stacked()`` differ only in how they resolve
        their per-panel lists and in per-panel policy: how many SNPs to
        label, the association panels' height and the floor on the figure
        height. Everything else is here.

        Args:
            config: Validated region, column, display, LD, and panel settings.
            association_inputs: Selected frames and resolved options, one per panel.
            threshold: Resolved p-value for the significance line, or None.
            label_top_n: SNPs to label per panel when the display config
                leaves it unset.
            association_height: Height-ratio units for each association panel.
            min_figure_height: Floor on the figure height in inches.
            recomb_lifter: Lifter for the plotter-loaded recombination maps,
                or None to use the registered chain.
        """
        region = config.region
        display = config.display.with_defaults(
            label_top_n=label_top_n, auto_genes=self._auto_genes
        )
        inputs = config.panels
        genes_df, exons_df = inputs.genes_df, inputs.exons_df
        recomb_df = inputs.recomb_df

        if genes_df is None and display.auto_genes:
            logger.debug(
                "auto_genes enabled, fetching genes for chr{}:{}-{}",
                region.chrom,
                region.start,
                region.end,
            )
            try:
                annotations = get_genes_for_build(
                    source_for(self.species, self.genome_build),
                    region.chrom,
                    region.start,
                    region.end,
                )
            except ReferenceAPIError as e:
                warnings.warn(
                    f"Gene track skipped for chr{region.chrom}:{region.start}-"
                    f"{region.end}; the gene source failed: {e}",
                    stacklevel=3,
                )
            else:
                if annotations.genes.empty:
                    logger.debug("No genes found in region")
                else:
                    genes_df = annotations.genes
                    if exons_df is None:
                        exons_df = annotations.exons

        finemap = (
            FinemappingPanel.from_frame(
                inputs.finemapping.data,
                region,
                inputs.finemapping.cs_col,
                chrom_col=inputs.finemapping.chrom_col,
            )
            if inputs.finemapping is not None
            else None
        )
        eqtl = (
            EqtlPanel.from_frame(
                inputs.eqtl.data,
                region,
                inputs.eqtl.gene,
                inputs.eqtl.threshold,
                chrom_col=inputs.eqtl.chrom_col,
            )
            if inputs.eqtl is not None
            else None
        )
        genes = (
            GenePanel.from_genes(genes_df, region, exons_df)
            if genes_df is not None
            else None
        )

        if display.show_recombination and recomb_df is None:
            recomb = self._get_recomb_for_region(
                region.chrom, region.start, region.end, recomb_lifter
            )
            if recomb.status is RecombStatus.OK:
                recomb_df = recomb.frame
            else:
                warnings.warn(
                    f"Recombination overlay skipped; {recomb.detail}",
                    stacklevel=3,
                )

        association: List[AssociationPanel] = []
        for index, request in enumerate(association_inputs):
            columns, ld = request.columns, request.ld
            df, ld_col = enrich_with_ld(
                request.data,
                reference_file=ld.ld_reference_file,
                lead_index=request.lead_index,
                ld_col=ld.ld_col,
                rs_col=request.rs_col,
                start=region.start,
                end=region.end,
                plink_path=self.plink_path,
                species=self.species,
                context=f"panel {index + 1}",
            )
            association.append(
                AssociationPanel(
                    data=df,
                    region=region,
                    height=association_height,
                    columns=columns,
                    display=display,
                    genomewide_threshold=threshold,
                    ld_col=ld_col,
                    hover=hover_for_association(columns, request.rs_col, ld_col),
                    lead_index=request.lead_index,
                    recomb_df=recomb_df if index == 0 else None,
                    panel_label=request.label,
                    add_ld_legend=(index == 0),
                )
            )

        heatmap = (
            HeatmapPanel.from_matrix(
                inputs.ld_heatmap.matrix,
                inputs.ld_heatmap.snp_ids,
                source=association[0],
                region=region,
                height=association_height * inputs.ld_heatmap.height,
                metric=inputs.ld_heatmap.metric,
            )
            if inputs.ld_heatmap is not None
            else None
        )

        panels: List[RegionalPanel] = [
            *association,
            *(panel for panel in (finemap, eqtl, genes, heatmap) if panel is not None),
        ]
        height = max(min_figure_height, sum(panel.height for panel in panels))
        logger.debug(
            "Creating regional plot with {} panels for chr{}:{}-{}",
            len(panels),
            region.chrom,
            region.start,
            region.end,
        )
        return render_figure(
            self._backend,
            FigurePlan(
                panels=panels,
                figsize=(display.figsize[0], height),
                height_ratios=[panel.height for panel in panels],
                xlabel=f"Chromosome {region.chrom} (Mb)",
                mb_xaxis=True,
                hspace=0.1,
            ),
        )

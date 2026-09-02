"""Main LocusZoomPlotter class for regional association plots.

Orchestrates all components (LD coloring, gene track, recombination overlay,
SNP labels) into a unified plotting interface.

Supports multiple backends:
- matplotlib (default): Static publication-quality plots
- plotly: Interactive HTML with hover tooltips
- bokeh: Interactive HTML for dashboards
"""

import warnings
from typing import Any, List, Optional, Union

import pandas as pd

from ._data import prepare_pvalue_data
from ._figure import FigurePlan, render_figure
from ._ld_plotting import enrich_with_ld
from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from ._regional_panels import (
    AssociationPanel,
    EqtlPanel,
    FinemappingPanel,
    GenePanel,
    HeatmapPanel,
    RegionalPanel,
    hover_for_association,
)
from .backends import BackendType, get_backend
from .config import (
    ColumnConfig,
    DisplayConfig,
    LDConfig,
    PanelInputs,
    PlotConfig,
    RegionConfig,
    StackedPlotConfig,
    resolve_deprecated_columns,
)
from .exceptions import ReferenceAPIError
from .ld import find_plink
from .logging import enable_logging, logger
from .recombination import RecombResult, RecombStatus, recomb_for_region
from .reference_genes import get_genes_for_build, source_for
from .schemas import validate_genes_df, validate_gwas_df
from .species import Species, resolve_species
from .utils import filter_by_region


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

    def _get_recomb_for_region(self, chrom: int, start: int, end: int) -> RecombResult:
        """Get a region's recombination rates, or the reason there are none.

        Caches per region and build. The caller renders the outcome; this does
        not warn, so a region asked for twice is reported once.
        """
        cache_key = (chrom, start, end, self.genome_build)
        if cache_key not in self._recomb_cache:
            self._recomb_cache[cache_key] = recomb_for_region(
                chrom=chrom,
                start=start,
                end=end,
                species=self.species,
                data_dir=self.recomb_data_dir,
                genome_build=self.genome_build,
            )
        return self._recomb_cache[cache_key]

    def plot(
        self,
        gwas_df: pd.DataFrame,
        *,
        chrom: Union[int, str],
        start: int,
        end: int,
        columns: ColumnConfig = ColumnConfig(),
        display: DisplayConfig = DisplayConfig(),
        ld: LDConfig = LDConfig(),
        panels: PanelInputs = PanelInputs(),
        significance_threshold: ThresholdArg = UNSET,
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

        Returns:
            Backend-specific figure object (``matplotlib.figure.Figure``,
            ``plotly.graph_objects.Figure``, or ``bokeh.layouts.Column``).

        Raises:
            ValueError: On an invalid region or a contradictory config
                (raised by :class:`PlotConfig` as a ``ValidationError``), or
                a missing required GWAS column.
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
        config = PlotConfig(
            region=RegionConfig(chrom=chrom, start=start, end=end),
            columns=columns,
            display=display,
            ld=ld,
            panels=panels,
        )
        return self._render_regional(
            config,
            [gwas_df],
            leads=None if ld.lead_pos is None else [ld.lead_pos],
            reference_files=[ld.ld_reference_file],
            panel_labels=None,
            threshold=resolve_threshold(
                significance_threshold, self.genomewide_threshold
            ),
            label_top_n=5,
            association_height=display.figsize[1] * 0.6,
            min_figure_height=0.0,
        )

    def plot_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
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
                the strongest in-region p-value when omitted. Required with
                a broadcast ``ld.ld_reference_file``.
            panel_labels: One label per panel, or None for none.
            ld_reference_files: One PLINK fileset per panel, replacing the
                broadcast ``ld.ld_reference_file``.
            significance_threshold: As on :meth:`plot`.

        Raises:
            ValueError: If ``gwas_dfs`` is empty or a per-panel list has a
                different length.

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
        if not gwas_dfs:
            raise ValueError("At least one GWAS DataFrame required")
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
        files = ld_reference_files or [ld.ld_reference_file] * len(gwas_dfs)
        return self._render_regional(
            config,
            gwas_dfs,
            leads=lead_positions,
            reference_files=files,
            panel_labels=panel_labels,
            threshold=resolve_threshold(
                significance_threshold, self.genomewide_threshold
            ),
            label_top_n=3,
            association_height=2.5,
            min_figure_height=display.figsize[1],
        )

    def _render_regional(
        self,
        config: PlotConfig,
        gwas_dfs: List[pd.DataFrame],
        *,
        leads: Optional[List[int]],
        reference_files: List[Optional[str]],
        panel_labels: Optional[List[str]],
        threshold: Optional[float],
        label_top_n: int,
        association_height: float,
        min_figure_height: float,
    ) -> Any:
        """Build the panel plan for one regional figure and render it.

        ``plot()`` and ``plot_stacked()`` differ only in how they resolve
        their per-panel lists and in per-panel policy: how many SNPs to
        label, the association panels' height and the floor on the figure
        height. Everything else is here.

        Args:
            config: Validated region, column, display, LD, and panel settings.
            gwas_dfs: One frame per association panel.
            leads: Lead position per panel, parallel to ``gwas_dfs``, or
                None to take each panel's strongest in-region signal.
            reference_files: PLINK fileset per panel, parallel to ``gwas_dfs``.
            panel_labels: Label per panel, or None for none.
            threshold: Resolved p-value for the significance line, or None.
            label_top_n: SNPs to label per panel when the display config
                leaves it unset.
            association_height: Height-ratio units for each association panel.
            min_figure_height: Floor on the figure height in inches.
        """
        region = config.region
        columns = resolve_deprecated_columns(gwas_dfs[0], config.columns)
        display = config.display.with_defaults(
            label_top_n=label_top_n, auto_genes=self._auto_genes
        )
        inputs = config.panels
        genes_df, exons_df = inputs.genes_df, inputs.exons_df
        recomb_df = inputs.recomb_df
        for gwas_df in gwas_dfs:
            validate_gwas_df(gwas_df, pos_col=columns.pos_col, p_col=columns.p_col)

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
        if genes_df is not None:
            validate_genes_df(genes_df)

        finemap = (
            FinemappingPanel.from_frame(
                inputs.finemapping_df, region, inputs.finemapping_cs_col
            )
            if inputs.finemapping_df is not None
            else None
        )
        eqtl = (
            EqtlPanel.from_frame(
                inputs.eqtl_df, region, inputs.eqtl_gene, inputs.eqtl_threshold
            )
            if inputs.eqtl_df is not None
            else None
        )
        genes = (
            GenePanel.from_genes(genes_df, region, exons_df)
            if genes_df is not None
            else None
        )

        if display.show_recombination and recomb_df is None:
            recomb = self._get_recomb_for_region(region.chrom, region.start, region.end)
            if recomb.status is RecombStatus.OK:
                recomb_df = recomb.frame
            else:
                warnings.warn(
                    f"Recombination overlay skipped; {recomb.detail}",
                    stacklevel=3,
                )

        association: List[AssociationPanel] = []
        for index, (gwas_df, reference_file) in enumerate(
            zip(gwas_dfs, reference_files)
        ):
            prepared = prepare_pvalue_data(gwas_df, columns.p_col)
            lead_pos = (
                leads[index]
                if leads is not None
                else _strongest_position(prepared, region, columns)
            )
            df, ld_col = enrich_with_ld(
                prepared,
                reference_file=reference_file,
                lead_pos=lead_pos,
                ld_col=config.ld.ld_col,
                pos_col=columns.pos_col,
                rs_col=columns.rs_col,
                start=region.start,
                end=region.end,
                plink_path=self.plink_path,
                species=self.species,
                context=f"panel {index + 1}",
            )
            if ld_col is not None and ld_col not in df.columns:
                ld_col = None
            association.append(
                AssociationPanel(
                    data=df,
                    region=region,
                    height=association_height,
                    columns=columns,
                    display=display,
                    genomewide_threshold=threshold,
                    ld_col=ld_col,
                    hover=hover_for_association(df, columns, ld_col),
                    lead_pos=lead_pos,
                    recomb_df=recomb_df if index == 0 else None,
                    panel_label=panel_labels[index] if panel_labels else None,
                    add_ld_legend=(index == 0),
                )
            )

        heatmap = (
            HeatmapPanel.from_matrix(
                inputs.ld_heatmap_df,
                inputs.ld_heatmap_snp_ids,
                source=association[0],
                region=region,
                height=association_height * inputs.ld_heatmap_height,
                metric=inputs.ld_heatmap_metric,
            )
            if inputs.ld_heatmap_df is not None
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


def _strongest_position(
    df: pd.DataFrame, region: RegionConfig, columns: ColumnConfig
) -> Optional[int]:
    """Return the position of the strongest in-region signal in a prepared frame.

    ``df`` has been through ``prepare_pvalue_data``, so the p-value domain rule
    has already been applied and ``neglog10p`` is finite. Filters on the
    canonical chromosome column when the frame carries one, so a whole-genome
    frame cannot anchor the lead to another chromosome.
    """
    region_df = filter_by_region(
        df,
        region=(region.chrom, region.start, region.end),
        pos_col=columns.pos_col,
    )
    if region_df.empty:
        logger.warning("No valid p-values in region, cannot determine lead SNP")
        return None
    return int(region_df.loc[region_df["neglog10p"].idxmax(), columns.pos_col])

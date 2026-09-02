"""Main LocusZoomPlotter class for regional association plots.

Orchestrates all components (LD coloring, gene track, recombination overlay,
SNP labels) into a unified plotting interface.

Supports multiple backends:
- matplotlib (default): Static publication-quality plots
- plotly: Interactive HTML with hover tooltips
- bokeh: Interactive HTML for dashboards
"""

import warnings
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ._data import prepare_pvalue_data
from ._ld_plotting import enrich_with_ld
from ._plotter_utils import DEFAULT_GENOMEWIDE_THRESHOLD
from ._regional import (
    AssociationPanel,
    EqtlPanel,
    FinemappingPanel,
    GenePanel,
    HeatmapPanel,
    RegionalFigurePlan,
    RegionalPanel,
    RegionalPlotComposer,
)
from .backends import BackendType, get_backend
from .config import ColumnConfig, PlotConfig, RegionConfig, StackedPlotConfig
from .exceptions import ReferenceAPIError
from .ld import find_plink
from .logging import enable_logging, logger
from .recombination import (
    ensure_recomb_maps,
    get_recombination_rate_for_region,
)
from .reference_genes import get_genes_for_build, source_for
from .utils import filter_by_region
from .validation import validate_genes_df, validate_gwas_df

# Precomputed significance line value (used for plotting)
DEFAULT_GENOMEWIDE_LINE = -np.log10(DEFAULT_GENOMEWIDE_THRESHOLD)


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
        species: Species name ('canine', 'feline', or None for custom).
            Canine has built-in recombination maps.
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
        >>> fig = plotter.plot(
        ...     gwas_df,
        ...     chrom=1,
        ...     start=1000000,
        ...     end=2000000,
        ...     lead_pos=1500000,
        ... )
        >>> fig.savefig("regional_plot.png", dpi=150)  # matplotlib
        >>> # or fig.save("plot.html")  # plotly/bokeh
    """

    def __init__(
        self,
        species: str = "canine",
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

        self.species = species
        self.genome_build = (
            genome_build if genome_build else self._default_build(species)
        )
        self._backend = get_backend(backend)
        self.plink_path = plink_path or find_plink()
        self.recomb_data_dir = recomb_data_dir
        self.genomewide_threshold = genomewide_threshold
        self._genomewide_line = -np.log10(genomewide_threshold)
        self._regional_composer = RegionalPlotComposer(
            self._backend,
            self._genomewide_line,
        )
        self._auto_genes = auto_genes
        self._recomb_cache = {}

    @staticmethod
    def _default_build(species: str) -> Optional[str]:
        """Get default genome build for species."""
        builds = {"canine": "canfam3.1", "feline": "felCat9"}
        return builds.get(species)

    def _ensure_recomb_maps(self) -> Optional[Path]:
        """Ensure recombination maps are available."""
        return ensure_recomb_maps(species=self.species, data_dir=self.recomb_data_dir)

    def _get_recomb_for_region(
        self, chrom: int, start: int, end: int
    ) -> Optional[pd.DataFrame]:
        """Get recombination rate data for a region, with caching."""
        cache_key = (chrom, start, end, self.genome_build)
        if cache_key in self._recomb_cache:
            return self._recomb_cache[cache_key]

        recomb_dir = self._ensure_recomb_maps()
        if recomb_dir is None:
            return None

        try:
            recomb_df = get_recombination_rate_for_region(
                chrom=chrom,
                start=start,
                end=end,
                species=self.species,
                data_dir=str(recomb_dir),
                genome_build=self.genome_build,
            )
            self._recomb_cache[cache_key] = recomb_df
            return recomb_df
        except FileNotFoundError as e:
            logger.warning(
                f"Recombination data file not found: {e}. "
                f"Recombination overlay will be skipped."
            )
            return None
        except ImportError as e:
            if "pyliftover" in str(e):
                import warnings

                warnings.warn(str(e), stacklevel=2)
                return None
            raise

    def plot(
        self,
        gwas_df: pd.DataFrame,
        *,
        chrom: Union[int, str],
        start: int,
        end: int,
        pos_col: str = "ps",
        p_col: str = "p_wald",
        rs_col: str = "rs",
        snp_labels: bool = True,
        label_top_n: int = 5,
        show_recombination: bool = True,
        figsize: Tuple[float, float] = (12.0, 8.0),
        lead_pos: Optional[int] = None,
        ld_reference_file: Optional[str] = None,
        ld_col: Optional[str] = None,
        genes_df: Optional[pd.DataFrame] = None,
        exons_df: Optional[pd.DataFrame] = None,
        recomb_df: Optional[pd.DataFrame] = None,
        eqtl_df: Optional[pd.DataFrame] = None,
        eqtl_gene: Optional[str] = None,
        eqtl_threshold: float = 1e-5,
        finemapping_df: Optional[pd.DataFrame] = None,
        finemapping_cs_col: Optional[str] = "cs",
        ld_heatmap_df: Optional[pd.DataFrame] = None,
        ld_heatmap_snp_ids: Optional[List[str]] = None,
        ld_heatmap_height: float = 0.25,
        ld_heatmap_metric: str = "r2",
    ) -> Any:
        """Create a regional association plot for a single locus.

        Plots ``-log10(p)`` against genomic position for the specified region,
        optionally overlaid with LD colouring, recombination rate, SNP labels,
        a gene track, and fine-mapping, eQTL, or LD-heatmap panels beneath.

        Coordinates are 1-based genomic positions (so ``lead_pos=0`` is
        rejected upstream by the Pydantic validator; any position ``>= 1`` is
        honoured, including small sentinel values).

        Args:
            gwas_df: GWAS summary statistics. Must contain the columns named
                by ``pos_col`` and ``p_col``; needs ``rs_col`` only when LD is
                computed from ``ld_reference_file``.
            chrom: Chromosome of the region.
            start: Region start position (bp, inclusive).
            end: Region end position (bp, inclusive).
            lead_pos: Genomic position of the lead SNP (``>= 1``). Required
                when ``ld_reference_file`` is supplied and ``ld_col`` is not.
            eqtl_df: eQTL results with ``pos`` and ``p_value`` columns; adds
                an eQTL panel. ``eqtl_gene`` filters it to one gene and
                requires a ``gene`` column; ``eqtl_threshold`` places its
                significance line.
            finemapping_df: Fine-mapping results with ``pos`` and ``pip``
                columns; adds a PIP panel coloured by ``finemapping_cs_col``.

        Returns:
            Backend-specific figure object (``matplotlib.figure.Figure``,
            ``plotly.graph_objects.Figure``, or ``bokeh.layouts.Column``).
            See keyword args for the full option surface.

        Raises:
            ValueError: On invalid kwargs (caught by
                :class:`PlotConfig`) or missing required GWAS columns.
            pylocuszoom.exceptions.PlinkError: When PLINK itself fails
                (timeout, non-zero exit, corrupt ``.bed``, missing output).
                The specific "empty LD output" case — singleton lead SNP with
                no neighbours in the window — is downgraded to a warning and
                the plot is drawn without LD colouring.
        """
        config = PlotConfig.from_kwargs(
            chrom=chrom,
            start=start,
            end=end,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
            snp_labels=snp_labels,
            label_top_n=label_top_n,
            show_recombination=show_recombination,
            figsize=figsize,
            lead_pos=lead_pos,
            ld_reference_file=ld_reference_file,
            ld_col=ld_col,
        )
        return self._render_regional(
            config,
            [gwas_df],
            leads=[config.ld.lead_pos],
            reference_files=[config.ld.ld_reference_file],
            panel_labels=None,
            association_height=config.display.figsize[1] * 0.6,
            min_figure_height=0.0,
            auto_genes=self._auto_genes,
            genes_df=genes_df,
            exons_df=exons_df,
            recomb_df=recomb_df,
            eqtl_df=eqtl_df,
            eqtl_gene=eqtl_gene,
            eqtl_threshold=eqtl_threshold,
            finemapping_df=finemapping_df,
            finemapping_cs_col=finemapping_cs_col,
            ld_heatmap_df=ld_heatmap_df,
            ld_heatmap_snp_ids=ld_heatmap_snp_ids,
            ld_heatmap_height=ld_heatmap_height,
            ld_heatmap_metric=ld_heatmap_metric,
        )

    def plot_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        *,
        chrom: Union[int, str],
        start: int,
        end: int,
        pos_col: str = "ps",
        p_col: str = "p_wald",
        rs_col: str = "rs",
        snp_labels: bool = True,
        label_top_n: int = 3,
        show_recombination: bool = True,
        figsize: Tuple[float, float] = (12.0, 8.0),
        ld_reference_file: Optional[str] = None,
        ld_col: Optional[str] = None,
        lead_positions: Optional[List[int]] = None,
        panel_labels: Optional[List[str]] = None,
        ld_reference_files: Optional[List[str]] = None,
        auto_genes: Optional[bool] = None,
        genes_df: Optional[pd.DataFrame] = None,
        exons_df: Optional[pd.DataFrame] = None,
        eqtl_df: Optional[pd.DataFrame] = None,
        eqtl_gene: Optional[str] = None,
        eqtl_threshold: float = 1e-5,
        finemapping_df: Optional[pd.DataFrame] = None,
        finemapping_cs_col: Optional[str] = "cs",
        recomb_df: Optional[pd.DataFrame] = None,
        ld_heatmap_df: Optional[pd.DataFrame] = None,
        ld_heatmap_snp_ids: Optional[List[str]] = None,
        ld_heatmap_height: float = 0.25,
        ld_heatmap_metric: str = "r2",
    ) -> Any:
        """Create stacked regional association plots for multiple GWAS.

        Each frame in ``gwas_dfs`` becomes one association panel; optional
        fine-mapping, eQTL, gene-track, and LD-heatmap panels follow beneath.
        ``lead_positions`` names one lead per panel and is auto-detected as
        the strongest in-region p-value when omitted.

        Args:
            gwas_dfs: One GWAS summary-statistics frame per panel.
            auto_genes: Fetch the gene track when ``genes_df`` is not given.
                ``None`` inherits the plotter's constructor setting.

        Raises:
            ValueError: If ``gwas_dfs`` is empty or a per-panel list
                (``lead_positions``, ``panel_labels``,
                ``ld_reference_files``) has a different length.
        """
        config = StackedPlotConfig.from_kwargs(
            chrom=chrom,
            start=start,
            end=end,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
            snp_labels=snp_labels,
            label_top_n=label_top_n,
            show_recombination=show_recombination,
            figsize=figsize,
            ld_reference_file=ld_reference_file,
            ld_col=ld_col,
            lead_positions=lead_positions,
            panel_labels=panel_labels,
            ld_reference_files=ld_reference_files,
        )
        n_gwas = len(gwas_dfs)
        if n_gwas == 0:
            raise ValueError("At least one GWAS DataFrame required")

        if config.lead_positions is not None and len(config.lead_positions) != n_gwas:
            raise ValueError(
                f"lead_positions length ({len(config.lead_positions)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )
        if config.panel_labels is not None and len(config.panel_labels) != n_gwas:
            raise ValueError(
                f"panel_labels length ({len(config.panel_labels)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )
        if (
            config.ld_reference_files is not None
            and len(config.ld_reference_files) != n_gwas
        ):
            raise ValueError(
                "ld_reference_files length "
                f"({len(config.ld_reference_files)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )

        leads = (
            list(config.lead_positions)
            if config.lead_positions is not None
            else _auto_lead_positions(gwas_dfs, config.region, config.columns)
        )
        files = config.ld_reference_files or [config.ld.ld_reference_file] * n_gwas
        return self._render_regional(
            config,
            gwas_dfs,
            leads=leads,
            reference_files=files,
            panel_labels=config.panel_labels,
            association_height=2.5,
            min_figure_height=config.display.figsize[1],
            auto_genes=self._auto_genes if auto_genes is None else auto_genes,
            genes_df=genes_df,
            exons_df=exons_df,
            recomb_df=recomb_df,
            eqtl_df=eqtl_df,
            eqtl_gene=eqtl_gene,
            eqtl_threshold=eqtl_threshold,
            finemapping_df=finemapping_df,
            finemapping_cs_col=finemapping_cs_col,
            ld_heatmap_df=ld_heatmap_df,
            ld_heatmap_snp_ids=ld_heatmap_snp_ids,
            ld_heatmap_height=ld_heatmap_height,
            ld_heatmap_metric=ld_heatmap_metric,
        )

    def _render_regional(
        self,
        config: PlotConfig,
        gwas_dfs: List[pd.DataFrame],
        *,
        leads: List[Optional[int]],
        reference_files: List[Optional[str]],
        panel_labels: Optional[List[str]],
        association_height: float,
        min_figure_height: float,
        auto_genes: bool,
        genes_df: Optional[pd.DataFrame],
        exons_df: Optional[pd.DataFrame],
        recomb_df: Optional[pd.DataFrame],
        eqtl_df: Optional[pd.DataFrame],
        eqtl_gene: Optional[str],
        eqtl_threshold: float,
        finemapping_df: Optional[pd.DataFrame],
        finemapping_cs_col: Optional[str],
        ld_heatmap_df: Optional[pd.DataFrame],
        ld_heatmap_snp_ids: Optional[List[str]],
        ld_heatmap_height: float,
        ld_heatmap_metric: str,
    ) -> Any:
        """Build the panel plan for one regional figure and render it.

        ``plot()`` and ``plot_stacked()`` differ only in how they resolve
        their per-panel lists and in height policy: the association panels'
        height and the floor on the figure height. Everything else is here.

        Args:
            config: Validated region, column, display, and LD settings.
            gwas_dfs: One frame per association panel.
            leads: Lead position per panel, parallel to ``gwas_dfs``.
            reference_files: PLINK fileset per panel, parallel to ``gwas_dfs``.
            panel_labels: Label per panel, or None for none.
            association_height: Height-ratio units for each association panel.
            min_figure_height: Floor on the figure height in inches.
            auto_genes: Fetch the gene track when ``genes_df`` is None.
        """
        region, columns, display = config.region, config.columns, config.display
        if ld_heatmap_df is not None and ld_heatmap_snp_ids is None:
            raise ValueError(
                "ld_heatmap_snp_ids is required when ld_heatmap_df is provided"
            )
        for gwas_df in gwas_dfs:
            validate_gwas_df(gwas_df, pos_col=columns.pos_col, p_col=columns.p_col)

        if genes_df is None and auto_genes:
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
            FinemappingPanel.from_frame(finemapping_df, region, finemapping_cs_col)
            if finemapping_df is not None
            else None
        )
        eqtl = (
            EqtlPanel.from_frame(eqtl_df, region, eqtl_gene, eqtl_threshold)
            if eqtl_df is not None
            else None
        )
        genes = (
            GenePanel.from_genes(genes_df, region, exons_df)
            if genes_df is not None
            else None
        )

        if display.show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(
                region.chrom, region.start, region.end
            )

        association: List[AssociationPanel] = []
        for index, (gwas_df, lead_pos, reference_file) in enumerate(
            zip(gwas_dfs, leads, reference_files)
        ):
            df, ld_col = enrich_with_ld(
                prepare_pvalue_data(gwas_df, columns.p_col),
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
            association.append(
                AssociationPanel(
                    data=df,
                    height=association_height,
                    columns=columns,
                    display=display,
                    ld_col=ld_col,
                    lead_pos=lead_pos,
                    recomb_df=recomb_df if index == 0 else None,
                    panel_label=panel_labels[index] if panel_labels else None,
                    add_ld_legend=(index == 0),
                )
            )

        heatmap = (
            HeatmapPanel.from_matrix(
                ld_heatmap_df,
                ld_heatmap_snp_ids,
                source=association[0],
                region=region,
                height=association_height * ld_heatmap_height,
                metric=ld_heatmap_metric,
            )
            if ld_heatmap_df is not None
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
        return self._regional_composer.render(
            RegionalFigurePlan(
                chrom=region.chrom,
                start=region.start,
                end=region.end,
                panels=panels,
                figsize=(display.figsize[0], height),
            )
        )


def _auto_lead_positions(
    gwas_dfs: List[pd.DataFrame], region: RegionConfig, columns: ColumnConfig
) -> List[Optional[int]]:
    """Pick each panel's lead as its strongest valid in-region p-value.

    Filters on a ``chrom`` or ``chr`` column when one exists, so a
    whole-genome frame cannot anchor the lead to another chromosome.
    Returns None for a panel with no valid in-region p-value.
    """
    leads: List[Optional[int]] = []
    for df in gwas_dfs:
        chrom_col = next((c for c in ("chrom", "chr") if c in df.columns), "chrom")
        region_df = filter_by_region(
            df,
            region=(region.chrom, region.start, region.end),
            chrom_col=chrom_col,
            pos_col=columns.pos_col,
        )
        if region_df.empty:
            leads.append(None)
            continue
        valid_p = region_df[columns.p_col].dropna()
        valid_p = valid_p[(valid_p >= 0) & (valid_p <= 1)]
        if valid_p.empty:
            logger.warning("No valid p-values in region, cannot determine lead SNP")
            leads.append(None)
        else:
            leads.append(int(region_df.loc[valid_p.idxmin(), columns.pos_col]))
    return leads

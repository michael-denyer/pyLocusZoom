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
from .config import PlotConfig, StackedPlotConfig
from .eqtl import validate_eqtl_df
from .exceptions import ReferenceAPIError
from .ld import find_plink
from .logging import enable_logging, logger
from .recombination import (
    ensure_recomb_maps,
    get_recombination_rate_for_region,
)
from .reference_genes import get_genes_for_build
from .utils import filter_by_region, validate_genes_df, validate_gwas_df

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
        ld_heatmap_df: Optional[pd.DataFrame] = None,
        ld_heatmap_snp_ids: Optional[List[str]] = None,
        ld_heatmap_height: float = 0.25,
        ld_heatmap_metric: str = "r2",
    ) -> Any:
        """Create a regional association plot for a single locus.

        Plots ``-log10(p)`` against genomic position for the specified region,
        optionally overlaid with LD colouring, recombination rate, SNP labels,
        a gene track, and/or a fine-mapping / LD-heatmap side panel.

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
        region = config.region
        columns = config.columns
        display = config.display
        ld_config = config.ld
        validate_gwas_df(
            gwas_df,
            pos_col=columns.pos_col,
            p_col=columns.p_col,
        )

        if ld_heatmap_df is not None and ld_heatmap_snp_ids is None:
            raise ValueError(
                "ld_heatmap_snp_ids is required when ld_heatmap_df is provided"
            )

        if genes_df is None and self._auto_genes:
            logger.debug(
                "auto_genes enabled, fetching genes for chr{}:{}-{}",
                region.chrom,
                region.start,
                region.end,
            )
            try:
                genes_df = get_genes_for_build(
                    species=self.species,
                    chrom=region.chrom,
                    start=region.start,
                    end=region.end,
                    genome_build=self.genome_build,
                    raise_on_error=True,
                )
            except ReferenceAPIError as e:
                warnings.warn(
                    f"Gene track skipped for chr{region.chrom}:{region.start}-"
                    f"{region.end}; the gene source failed: {e}",
                    stacklevel=2,
                )
                genes_df = None
            else:
                if genes_df.empty:
                    logger.debug("No genes found in region")
                    genes_df = None

        if genes_df is not None:
            validate_genes_df(genes_df)

        logger.debug(f"Creating plot for chr{region.chrom}:{region.start}-{region.end}")
        # Don't flip global matplotlib interactive mode here — it leaks into
        # the caller's notebook session and disables auto-display for all
        # subsequent plots from any library. Backends manage their own state.

        df = prepare_pvalue_data(gwas_df, columns.p_col)

        df, resolved_ld_col = enrich_with_ld(
            df,
            reference_file=ld_config.ld_reference_file,
            lead_pos=ld_config.lead_pos,
            ld_col=ld_config.ld_col,
            pos_col=columns.pos_col,
            rs_col=columns.rs_col,
            start=region.start,
            end=region.end,
            plink_path=self.plink_path,
            species=self.species,
        )

        if display.show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(
                region.chrom, region.start, region.end
            )

        association_height = display.figsize[1] * 0.6
        association = AssociationPanel(
            data=df,
            height=association_height,
            columns=columns,
            display=display,
            ld_col=resolved_ld_col,
            lead_pos=ld_config.lead_pos,
            recomb_df=recomb_df,
            add_ld_legend=True,
        )
        panels: List[RegionalPanel] = [association]

        if genes_df is not None:
            panels.append(GenePanel.from_genes(genes_df, region, exons_df))

        if ld_heatmap_df is not None:
            heatmap = HeatmapPanel.from_matrix(
                ld_heatmap_df,
                ld_heatmap_snp_ids,
                source=association,
                region=region,
                height=association_height * ld_heatmap_height,
                metric=ld_heatmap_metric,
            )
            if heatmap is not None:
                panels.append(heatmap)

        return self._regional_composer.render(
            RegionalFigurePlan(
                chrom=region.chrom,
                start=region.start,
                end=region.end,
                panels=panels,
                figsize=(display.figsize[0], sum(panel.height for panel in panels)),
            )
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
        """Create stacked regional association plots for multiple GWAS."""
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
        region = config.region
        columns = config.columns
        display = config.display
        ld_config = config.ld

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

        for df in gwas_dfs:
            validate_gwas_df(
                df,
                pos_col=columns.pos_col,
                p_col=columns.p_col,
            )
        if genes_df is not None:
            validate_genes_df(genes_df)
        if eqtl_df is not None:
            validate_eqtl_df(eqtl_df)

        if ld_heatmap_df is not None and ld_heatmap_snp_ids is None:
            raise ValueError(
                "ld_heatmap_snp_ids is required when ld_heatmap_df is provided"
            )

        resolved_lead_positions: List[Optional[int]] = (
            list(config.lead_positions) if config.lead_positions is not None else []
        )
        if config.lead_positions is None:
            for df in gwas_dfs:
                # Detect a chrom column under either common convention
                # ("chrom" or "chr") to prevent cross-chromosome lead
                # selection on whole-genome summary stats DataFrames.
                chrom_col = next(
                    (c for c in ("chrom", "chr") if c in df.columns), "chrom"
                )
                region_df = filter_by_region(
                    df,
                    region=(region.chrom, region.start, region.end),
                    chrom_col=chrom_col,
                    pos_col=columns.pos_col,
                )
                if not region_df.empty:
                    valid_p = region_df[columns.p_col].dropna()
                    valid_p = valid_p[(valid_p >= 0) & (valid_p <= 1)]
                    if valid_p.empty:
                        logger.warning(
                            "No valid p-values in region, cannot determine lead SNP"
                        )
                        resolved_lead_positions.append(None)
                    else:
                        lead_idx = valid_p.idxmin()
                        resolved_lead_positions.append(
                            int(region_df.loc[lead_idx, columns.pos_col])
                        )
                else:
                    resolved_lead_positions.append(None)

        resolved_reference_files = config.ld_reference_files
        if resolved_reference_files is None and ld_config.ld_reference_file is not None:
            resolved_reference_files = [ld_config.ld_reference_file] * n_gwas

        if display.show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(
                region.chrom, region.start, region.end
            )

        panel_height = 2.5
        association: List[AssociationPanel] = []
        for index, (gwas_df, lead_pos) in enumerate(
            zip(gwas_dfs, resolved_lead_positions)
        ):
            df = prepare_pvalue_data(gwas_df, columns.p_col)
            reference_file = (
                resolved_reference_files[index] if resolved_reference_files else None
            )
            df, panel_ld_col = enrich_with_ld(
                df,
                reference_file=reference_file,
                lead_pos=lead_pos,
                ld_col=ld_config.ld_col,
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
                    height=panel_height,
                    columns=columns,
                    display=display,
                    ld_col=panel_ld_col,
                    lead_pos=lead_pos,
                    recomb_df=recomb_df if index == 0 else None,
                    panel_label=(
                        config.panel_labels[index] if config.panel_labels else None
                    ),
                    add_ld_legend=(index == 0),
                )
            )
        panels: List[RegionalPanel] = list(association)

        if finemapping_df is not None:
            panels.append(
                FinemappingPanel.from_frame(finemapping_df, region, finemapping_cs_col)
            )

        if eqtl_df is not None:
            eqtl_data = eqtl_df.copy()

            eqtl_gene_filtered = False
            if eqtl_gene:
                if "gene" in eqtl_data.columns:
                    eqtl_data = eqtl_data[eqtl_data["gene"] == eqtl_gene]
                    eqtl_gene_filtered = True
                else:
                    logger.warning(
                        f"eqtl_gene='{eqtl_gene}' specified but eQTL data has no 'gene' column; "
                        "showing all eQTL data unfiltered"
                    )

            if "pos" in eqtl_data.columns:
                mask = (eqtl_data["pos"] >= region.start) & (
                    eqtl_data["pos"] <= region.end
                )
                if "chr" in eqtl_data.columns:
                    chrom_str = str(region.chrom).replace("chr", "")
                    eqtl_chrom = (
                        eqtl_data["chr"].astype(str).str.replace("chr", "", regex=False)
                    )
                    mask = mask & (eqtl_chrom == chrom_str)
                eqtl_data = eqtl_data[mask]

            if not eqtl_data.empty:
                eqtl_data = prepare_pvalue_data(eqtl_data, "p_value")
            panels.append(
                EqtlPanel(
                    data=eqtl_data,
                    height=2.0,
                    gene_filtered=eqtl_gene_filtered,
                    gene=eqtl_gene,
                    threshold=eqtl_threshold,
                )
            )

        if genes_df is not None:
            panels.append(GenePanel.from_genes(genes_df, region, exons_df))

        if ld_heatmap_df is not None:
            heatmap = HeatmapPanel.from_matrix(
                ld_heatmap_df,
                ld_heatmap_snp_ids,
                source=association[0],
                region=region,
                height=panel_height * ld_heatmap_height,
                metric=ld_heatmap_metric,
            )
            if heatmap is not None:
                panels.append(heatmap)

        total_height = max(
            display.figsize[1],
            sum(panel.height for panel in panels),
        )
        logger.debug(
            "Creating stacked plot with {} panels for chr{}:{}-{}",
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
                figsize=(display.figsize[0], total_height),
            )
        )

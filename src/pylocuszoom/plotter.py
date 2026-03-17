"""Main LocusZoomPlotter class for regional association plots.

Orchestrates all components (LD coloring, gene track, recombination overlay,
SNP labels) into a unified plotting interface.

Supports multiple backends:
- matplotlib (default): Static publication-quality plots
- plotly: Interactive HTML with hover tooltips
- bokeh: Interactive HTML for dashboards
"""

from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    calculate_gene_track_rows,
    transform_pvalues,
)
from .backends import BackendType, get_backend
from .backends.hover import HoverConfig, HoverDataBuilder
from .colors import (
    EQTL_NEGATIVE_BINS,
    EQTL_POSITIVE_BINS,
    LD_BINS,
    LD_HEATMAP_COLORS,
    LEAD_SNP_COLOR,
    LEAD_SNP_HIGHLIGHT_COLOR,
    get_credible_set_color,
    get_eqtl_color,
    get_ld_bin,
    get_ld_color_palette,
)
from .config import PlotConfig, StackedPlotConfig
from .ensembl import get_genes_for_region
from .eqtl import validate_eqtl_df
from .finemapping import (
    get_credible_sets,
    plot_finemapping,
    prepare_finemapping_for_plotting,
)
from .gene_track import plot_gene_track_generic
from .ld import calculate_ld, find_plink
from .logging import enable_logging, logger
from .recombination import (
    ensure_recomb_maps,
    get_recombination_rate_for_region,
)
from .utils import validate_genes_df, validate_gwas_df

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
        self._backend_name = backend
        self.plink_path = plink_path or find_plink()
        self.recomb_data_dir = recomb_data_dir
        self.genomewide_threshold = genomewide_threshold
        self._genomewide_line = -np.log10(genomewide_threshold)
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
        """Create a regional association plot."""
        PlotConfig.from_kwargs(
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
        validate_gwas_df(gwas_df, pos_col=pos_col, p_col=p_col)

        if ld_heatmap_df is not None and ld_heatmap_snp_ids is None:
            raise ValueError(
                "ld_heatmap_snp_ids is required when ld_heatmap_df is provided"
            )

        if genes_df is None and self._auto_genes:
            logger.debug(
                f"auto_genes enabled, fetching genes for chr{chrom}:{start}-{end}"
            )
            genes_df = get_genes_for_region(
                species=self.species,
                chrom=chrom,
                start=start,
                end=end,
            )
            if genes_df.empty:
                logger.debug("No genes found in region from Ensembl")
                genes_df = None

        if genes_df is not None:
            validate_genes_df(genes_df)

        logger.debug(f"Creating plot for chr{chrom}:{start}-{end}")
        plt.ioff()

        df = transform_pvalues(gwas_df, p_col)

        if ld_reference_file and lead_pos and ld_col is None:
            if rs_col not in df.columns:
                logger.warning(
                    f"Cannot calculate LD: column '{rs_col}' not found in GWAS data. "
                    f"Provide rs_col parameter or add SNP IDs to DataFrame."
                )
            else:
                lead_snp_row = df[df[pos_col] == lead_pos]
                if lead_snp_row.empty:
                    logger.warning(
                        f"Lead SNP at position {lead_pos} not found in GWAS data. "
                        f"LD coloring will be skipped."
                    )
                else:
                    lead_snp_id = lead_snp_row[rs_col].iloc[0]
                    logger.debug(f"Calculating LD for lead SNP {lead_snp_id}")
                    ld_df = calculate_ld(
                        bfile_path=ld_reference_file,
                        lead_snp=lead_snp_id,
                        window_kb=max((end - start) // 1000, 500),
                        plink_path=self.plink_path,
                        species=self.species,
                    )
                    if not ld_df.empty:
                        df = df.merge(
                            ld_df,
                            left_on=rs_col,
                            right_on="SNP",
                            how="left",
                            validate="many_to_one",
                        )
                        ld_col = "R2"

        if show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(chrom, start, end)

        heatmap_data = None
        if ld_heatmap_df is not None and ld_heatmap_snp_ids is not None:
            heatmap_data = self._transform_heatmap_to_genomic_coords(
                ld_matrix=ld_heatmap_df,
                snp_ids=ld_heatmap_snp_ids,
                gwas_df=df,
                start=start,
                end=end,
                rs_col=rs_col,
                pos_col=pos_col,
            )
            if heatmap_data is None:
                logger.warning(
                    "No SNPs from LD heatmap overlap with region - heatmap not rendered"
                )

        fig, ax, gene_ax, heatmap_ax = self._create_figure_with_heatmap(
            genes_df=genes_df,
            chrom=chrom,
            start=start,
            end=end,
            figsize=figsize,
            heatmap_data=heatmap_data,
            heatmap_height=ld_heatmap_height,
        )

        self._plot_association(ax, df, pos_col, ld_col, lead_pos, rs_col, p_col)

        self._backend.axhline(
            ax,
            y=self._genomewide_line,
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=0.65,
            zorder=1,
        )

        snp_label_texts: list = []
        if snp_labels and rs_col in df.columns and label_top_n > 0 and not df.empty:
            if self._backend.supports_snp_labels:
                snp_label_texts = self._backend.add_snp_labels(
                    ax,
                    df,
                    pos_col=pos_col,
                    neglog10p_col="neglog10p",
                    rs_col=rs_col,
                    label_top_n=label_top_n,
                    genes_df=genes_df,
                    chrom=chrom,
                    adjust=False,
                )

        if recomb_df is not None and not recomb_df.empty:
            if self._backend.supports_secondary_axis:
                self._add_recombination_overlay(ax, recomb_df, start, end)

        self._backend.set_ylabel(ax, r"$-\log_{10}$ P")
        self._backend.set_xlim(ax, start, end)
        has_recomb = recomb_df is not None and not recomb_df.empty
        if has_recomb and self._backend.supports_secondary_axis:
            self._backend.hide_spines(ax, ["top"])
        else:
            self._backend.hide_spines(ax, ["top", "right"])

        if ld_col is not None and ld_col in df.columns:
            self._backend.add_ld_legend(ax, LD_BINS, LEAD_SNP_COLOR)

        if genes_df is not None and gene_ax is not None:
            plot_gene_track_generic(
                gene_ax, self._backend, genes_df, chrom, start, end, exons_df
            )
            self._backend.hide_spines(gene_ax, ["top", "right", "left"])
            self._backend.format_xaxis_mb(gene_ax)
            if heatmap_ax is None:
                self._backend.set_xlabel(gene_ax, f"Chromosome {chrom} (Mb)")

        if heatmap_ax is not None and heatmap_data is not None:
            filtered_matrix, x_positions, filtered_snp_ids = heatmap_data
            lead_snp_id = None
            if lead_pos is not None and rs_col in df.columns:
                lead_row = df[df[pos_col] == lead_pos]
                if not lead_row.empty:
                    lead_snp_id = lead_row[rs_col].iloc[0]
            self._render_heatmap_panel(
                ax=heatmap_ax,
                fig=fig,
                ld_matrix=filtered_matrix,
                x_positions=x_positions,
                snp_ids=filtered_snp_ids,
                metric=ld_heatmap_metric,
                lead_snp_id=lead_snp_id,
                start=start,
                end=end,
            )
            self._backend.set_xlabel(heatmap_ax, f"Chromosome {chrom} (Mb)")
            self._backend.format_xaxis_mb(heatmap_ax)
        elif gene_ax is None and heatmap_ax is None:
            self._backend.set_xlabel(ax, f"Chromosome {chrom} (Mb)")

        self._backend.format_xaxis_mb(ax)
        self._backend.finalize_layout(fig, hspace=0.1)

        if snp_label_texts:
            self._backend.adjust_snp_labels(ax, snp_label_texts)

        return fig

    def _create_figure_with_heatmap(
        self,
        genes_df: Optional[pd.DataFrame],
        chrom: int,
        start: int,
        end: int,
        figsize: Tuple[float, float],
        heatmap_data: Optional[Tuple[pd.DataFrame, List[int], List[str]]],
        heatmap_height: float = 0.25,
    ) -> Tuple[Any, Any, Optional[Any], Optional[Any]]:
        """Create figure with optional gene track and heatmap panel."""
        assoc_height = figsize[1] * 0.6

        gene_track_height = 0.0
        if genes_df is not None:
            n_gene_rows = calculate_gene_track_rows(genes_df, chrom, start, end)
            base_gene_height = 1.0
            per_row_height = 0.5
            gene_track_height = base_gene_height + (n_gene_rows - 1) * per_row_height

        actual_heatmap_height = 0.0
        if heatmap_data is not None:
            actual_heatmap_height = assoc_height * heatmap_height

        n_panels = 1
        height_ratios = [assoc_height]

        if genes_df is not None:
            n_panels += 1
            height_ratios.append(gene_track_height)

        if heatmap_data is not None:
            n_panels += 1
            height_ratios.append(actual_heatmap_height)

        total_height = sum(height_ratios)

        fig, axes = self._backend.create_figure(
            n_panels=n_panels,
            height_ratios=height_ratios,
            figsize=(figsize[0], total_height),
            sharex=True,
        )

        assoc_ax = axes[0]
        gene_ax = None
        heatmap_ax = None

        panel_idx = 1
        if genes_df is not None:
            gene_ax = axes[panel_idx]
            panel_idx += 1
        if heatmap_data is not None:
            heatmap_ax = axes[panel_idx]

        return fig, assoc_ax, gene_ax, heatmap_ax

    def _transform_heatmap_to_genomic_coords(
        self,
        ld_matrix: pd.DataFrame,
        snp_ids: List[str],
        gwas_df: pd.DataFrame,
        start: int,
        end: int,
        rs_col: str,
        pos_col: str,
    ) -> Optional[Tuple[pd.DataFrame, List[int], List[str]]]:
        """Transform heatmap matrix to genomic coordinates."""
        if rs_col not in gwas_df.columns:
            logger.warning(
                f"Cannot map heatmap to genomic coords: column '{rs_col}' not in GWAS data"
            )
            return None

        snp_to_pos = dict(zip(gwas_df[rs_col], gwas_df[pos_col]))

        filtered_indices = []
        filtered_snp_ids = []
        x_positions = []

        for i, snp_id in enumerate(snp_ids):
            if snp_id in snp_to_pos:
                pos = snp_to_pos[snp_id]
                if start <= pos <= end:
                    filtered_indices.append(i)
                    filtered_snp_ids.append(snp_id)
                    x_positions.append(int(pos))

        if not filtered_indices:
            return None

        filtered_matrix = ld_matrix.iloc[filtered_indices, filtered_indices].copy()
        return filtered_matrix, x_positions, filtered_snp_ids

    def _render_heatmap_panel(
        self,
        ax: Any,
        fig: Any,
        ld_matrix: pd.DataFrame,
        x_positions: List[int],
        snp_ids: List[str],
        metric: str,
        lead_snp_id: Optional[str],
        start: int,
        end: int,
    ) -> None:
        """Render LD heatmap panel with genomic x-coordinates."""
        data = ld_matrix.values
        n_snps = len(snp_ids)

        if n_snps < 2:
            logger.debug("Skipping heatmap: fewer than 2 SNPs after filtering")
            return

        mappable = self._backend.add_heatmap(
            ax,
            data=data,
            x_coords=x_positions,
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
            mask_upper=True,
        )

        label = "R\u00b2" if metric == "r2" else "D'"
        self._backend.add_colorbar(ax, mappable, label=label)

        if lead_snp_id is not None and lead_snp_id in snp_ids:
            lead_idx = snp_ids.index(lead_snp_id)
            self._highlight_heatmap_snp(ax, fig, lead_idx, n_snps)

        self._backend.set_xlim(ax, start, end)
        self._backend.set_yticks(ax, [], [])
        self._backend.hide_spines(ax, ["top", "right", "left"])

    def _highlight_heatmap_snp(
        self, ax: Any, fig: Any, snp_idx: int, n_snps: int
    ) -> None:
        """Highlight a SNP's row/column in the heatmap."""
        self._backend.highlight_heatmap_snp(
            ax,
            fig,
            snp_idx,
            n_snps,
            color=LEAD_SNP_HIGHLIGHT_COLOR,
            linewidth=2,
        )

    def _plot_association(
        self,
        ax: Any,
        df: pd.DataFrame,
        pos_col: str,
        ld_col: Optional[str],
        lead_pos: Optional[int],
        rs_col: Optional[str] = None,
        p_col: Optional[str] = None,
    ) -> None:
        """Plot association scatter with LD coloring."""
        hover_config = HoverConfig(
            snp_col=rs_col if rs_col and rs_col in df.columns else None,
            pos_col=pos_col if pos_col in df.columns else None,
            p_col=p_col if p_col and p_col in df.columns else None,
            ld_col=ld_col if ld_col and ld_col in df.columns else None,
        )
        hover_builder = HoverDataBuilder(hover_config)

        if ld_col is not None and ld_col in df.columns:
            df["ld_bin"] = df[ld_col].apply(get_ld_bin)
            df = df.sort_values(ld_col, ascending=True, na_position="first")

            palette = get_ld_color_palette()
            for bin_label in df["ld_bin"].unique():
                bin_data = df[df["ld_bin"] == bin_label]
                self._backend.scatter(
                    ax,
                    bin_data[pos_col],
                    bin_data["neglog10p"],
                    colors=palette.get(bin_label, "#BEBEBE"),
                    sizes=60,
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                    hover_data=hover_builder.build_dataframe(bin_data),
                )
        else:
            self._backend.scatter(
                ax,
                df[pos_col],
                df["neglog10p"],
                colors="#BEBEBE",
                sizes=60,
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
                hover_data=hover_builder.build_dataframe(df),
            )

        if lead_pos is not None:
            lead_snp = df[df[pos_col] == lead_pos]
            if not lead_snp.empty:
                self._backend.scatter(
                    ax,
                    lead_snp[pos_col],
                    lead_snp["neglog10p"],
                    colors=LEAD_SNP_COLOR,
                    sizes=120,
                    marker="D",
                    edgecolor="black",
                    linewidth=1.5,
                    zorder=10,
                    hover_data=hover_builder.build_dataframe(lead_snp),
                )

    def _add_recombination_overlay(
        self,
        ax: Any,
        recomb_df: pd.DataFrame,
        start: int,
        end: int,
    ) -> None:
        """Add recombination overlay for all backends."""
        if not hasattr(self._backend, "add_recombination_overlay"):
            logger.warning(
                "Backend '{}' does not implement add_recombination_overlay, "
                "skipping recombination overlay",
                self._backend_name,
            )
            return
        self._backend.add_recombination_overlay(ax, recomb_df, start, end)

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
        StackedPlotConfig.from_kwargs(
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

        if lead_positions is not None and len(lead_positions) != n_gwas:
            raise ValueError(
                f"lead_positions length ({len(lead_positions)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )
        if panel_labels is not None and len(panel_labels) != n_gwas:
            raise ValueError(
                f"panel_labels length ({len(panel_labels)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )
        if ld_reference_files is not None and len(ld_reference_files) != n_gwas:
            raise ValueError(
                f"ld_reference_files length ({len(ld_reference_files)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )

        for i, df in enumerate(gwas_dfs):
            validate_gwas_df(df, pos_col=pos_col, p_col=p_col)
        if genes_df is not None:
            validate_genes_df(genes_df)
        if eqtl_df is not None:
            validate_eqtl_df(eqtl_df)

        if ld_heatmap_df is not None and ld_heatmap_snp_ids is None:
            raise ValueError(
                "ld_heatmap_snp_ids is required when ld_heatmap_df is provided"
            )

        if lead_positions is None:
            lead_positions = []
            for df in gwas_dfs:
                region_df = df[(df[pos_col] >= start) & (df[pos_col] <= end)]
                if not region_df.empty:
                    valid_p = region_df[p_col].dropna()
                    valid_p = valid_p[(valid_p >= 0) & (valid_p <= 1)]
                    if valid_p.empty:
                        logger.warning(
                            "No valid p-values in region, cannot determine lead SNP"
                        )
                        lead_positions.append(None)
                    else:
                        lead_idx = valid_p.idxmin()
                        lead_positions.append(int(region_df.loc[lead_idx, pos_col]))
                else:
                    lead_positions.append(None)

        if ld_reference_files is None and ld_reference_file is not None:
            ld_reference_files = [ld_reference_file] * n_gwas

        heatmap_data = None
        if ld_heatmap_df is not None and ld_heatmap_snp_ids is not None:
            first_gwas = transform_pvalues(gwas_dfs[0], p_col)
            heatmap_data = self._transform_heatmap_to_genomic_coords(
                ld_matrix=ld_heatmap_df,
                snp_ids=ld_heatmap_snp_ids,
                gwas_df=first_gwas,
                start=start,
                end=end,
                rs_col=rs_col,
                pos_col=pos_col,
            )
            if heatmap_data is None:
                logger.warning(
                    "No SNPs from LD heatmap overlap with region - heatmap not rendered"
                )

        panel_height = 2.5
        eqtl_height = 2.0 if eqtl_df is not None else 0
        finemapping_height = 1.5 if finemapping_df is not None else 0
        heatmap_height_inches = panel_height * ld_heatmap_height if heatmap_data else 0

        if genes_df is not None:
            n_gene_rows = calculate_gene_track_rows(genes_df, chrom, start, end)
            gene_track_height = 1.0 + (n_gene_rows - 1) * 0.5
        else:
            gene_track_height = 0

        n_panels = (
            n_gwas
            + (1 if finemapping_df is not None else 0)
            + (1 if eqtl_df is not None else 0)
            + (1 if genes_df is not None else 0)
            + (1 if heatmap_data is not None else 0)
        )
        height_ratios = [panel_height] * n_gwas
        if finemapping_df is not None:
            height_ratios.append(finemapping_height)
        if eqtl_df is not None:
            height_ratios.append(eqtl_height)
        if genes_df is not None:
            height_ratios.append(gene_track_height)
        if heatmap_data is not None:
            height_ratios.append(heatmap_height_inches)

        total_height = max(figsize[1], sum(height_ratios))
        actual_figsize = (figsize[0], total_height)

        logger.debug(
            f"Creating stacked plot with {n_panels} panels for chr{chrom}:{start}-{end}"
        )

        if show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(chrom, start, end)

        fig, axes = self._backend.create_figure(
            n_panels=n_panels,
            height_ratios=height_ratios,
            figsize=actual_figsize,
            sharex=True,
        )

        all_snp_label_texts: list[tuple] = []

        for i, (gwas_df, lead_pos) in enumerate(zip(gwas_dfs, lead_positions)):
            ax = axes[i]
            df = transform_pvalues(gwas_df, p_col)

            panel_ld_col = ld_col
            if ld_reference_files and ld_reference_files[i] and lead_pos and not ld_col:
                if rs_col not in df.columns:
                    logger.warning(
                        f"Cannot calculate LD for panel {i + 1}: column '{rs_col}' "
                        f"not found in GWAS data. "
                        f"Provide rs_col parameter or add SNP IDs to DataFrame."
                    )
                else:
                    lead_snp_row = df[df[pos_col] == lead_pos]
                    if not lead_snp_row.empty:
                        lead_snp_id = lead_snp_row[rs_col].iloc[0]
                        ld_df = calculate_ld(
                            bfile_path=ld_reference_files[i],
                            lead_snp=lead_snp_id,
                            window_kb=max((end - start) // 1000, 500),
                            plink_path=self.plink_path,
                            species=self.species,
                        )
                        if not ld_df.empty:
                            df = df.merge(
                                ld_df,
                                left_on=rs_col,
                                right_on="SNP",
                                how="left",
                                validate="many_to_one",
                            )
                            panel_ld_col = "R2"

            self._plot_association(
                ax, df, pos_col, panel_ld_col, lead_pos, rs_col, p_col
            )

            self._backend.axhline(
                ax,
                y=self._genomewide_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.65,
                zorder=1,
            )

            if snp_labels and rs_col in df.columns and label_top_n > 0 and not df.empty:
                if self._backend.supports_snp_labels:
                    texts = self._backend.add_snp_labels(
                        ax,
                        df,
                        pos_col=pos_col,
                        neglog10p_col="neglog10p",
                        rs_col=rs_col,
                        label_top_n=label_top_n,
                        genes_df=genes_df,
                        chrom=chrom,
                        adjust=False,
                    )
                    if texts:
                        all_snp_label_texts.append((ax, texts))

            if i == 0 and recomb_df is not None and not recomb_df.empty:
                if self._backend.supports_secondary_axis:
                    self._add_recombination_overlay(ax, recomb_df, start, end)

            self._backend.set_ylabel(ax, r"$-\log_{10}$ P")
            self._backend.set_xlim(ax, start, end)
            self._backend.hide_spines(ax, ["top", "right"])

            if panel_labels and i < len(panel_labels):
                self._backend.add_panel_label(ax, panel_labels[i])

            if i == 0 and panel_ld_col is not None and panel_ld_col in df.columns:
                self._backend.add_ld_legend(ax, LD_BINS, LEAD_SNP_COLOR)

        panel_idx = n_gwas

        if finemapping_df is not None:
            ax = axes[panel_idx]
            fm_data = prepare_finemapping_for_plotting(
                finemapping_df,
                pos_col="pos",
                pip_col="pip",
                chrom=chrom,
                start=start,
                end=end,
            )
            if not fm_data.empty:
                plot_finemapping(
                    self._backend,
                    ax,
                    fm_data,
                    pos_col="pos",
                    pip_col="pip",
                    cs_col=finemapping_cs_col,
                    show_credible_sets=True,
                    pip_threshold=0.01,
                )
                credible_sets = get_credible_sets(fm_data, finemapping_cs_col)
                if credible_sets:
                    self._backend.add_finemapping_legend(
                        ax, credible_sets, get_credible_set_color
                    )
            self._backend.set_ylabel(ax, "PIP")
            self._backend.set_ylim(ax, -0.05, 1.05)
            self._backend.hide_spines(ax, ["top", "right"])
            panel_idx += 1

        eqtl_panel_idx = panel_idx
        if eqtl_df is not None:
            ax = axes[eqtl_panel_idx]
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
                mask = (eqtl_data["pos"] >= start) & (eqtl_data["pos"] <= end)
                if "chr" in eqtl_data.columns:
                    chrom_str = str(chrom).replace("chr", "")
                    eqtl_chrom = (
                        eqtl_data["chr"].astype(str).str.replace("chr", "", regex=False)
                    )
                    mask = mask & (eqtl_chrom == chrom_str)
                eqtl_data = eqtl_data[mask]

            if not eqtl_data.empty:
                eqtl_data = transform_pvalues(eqtl_data, "p_value")

                eqtl_extra_cols = {}
                if "effect_size" in eqtl_data.columns:
                    eqtl_extra_cols["effect_size"] = "Effect"
                if "gene" in eqtl_data.columns:
                    eqtl_extra_cols["gene"] = "Gene"
                eqtl_hover_config = HoverConfig(
                    pos_col="pos" if "pos" in eqtl_data.columns else None,
                    p_col="p_value" if "p_value" in eqtl_data.columns else None,
                    extra_cols=eqtl_extra_cols,
                )
                eqtl_hover_builder = HoverDataBuilder(eqtl_hover_config)

                has_effect = "effect_size" in eqtl_data.columns

                if has_effect:
                    pos_effects = eqtl_data[eqtl_data["effect_size"] >= 0]
                    neg_effects = eqtl_data[eqtl_data["effect_size"] < 0]

                    if not pos_effects.empty:
                        pos_colors = pos_effects["effect_size"].apply(get_eqtl_color)
                        self._backend.scatter(
                            ax,
                            pos_effects["pos"],
                            pos_effects["neglog10p"],
                            colors=pos_colors.tolist(),
                            sizes=50,
                            marker="^",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                            hover_data=eqtl_hover_builder.build_dataframe(pos_effects),
                        )

                    if not neg_effects.empty:
                        neg_colors = neg_effects["effect_size"].apply(get_eqtl_color)
                        self._backend.scatter(
                            ax,
                            neg_effects["pos"],
                            neg_effects["neglog10p"],
                            colors=neg_colors.tolist(),
                            sizes=50,
                            marker="v",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                            hover_data=eqtl_hover_builder.build_dataframe(neg_effects),
                        )

                    self._backend.add_eqtl_legend(
                        ax, EQTL_POSITIVE_BINS, EQTL_NEGATIVE_BINS
                    )
                else:
                    label = f"eQTL ({eqtl_gene})" if eqtl_gene_filtered else "eQTL"
                    self._backend.scatter(
                        ax,
                        eqtl_data["pos"],
                        eqtl_data["neglog10p"],
                        colors="#FF6B6B",
                        sizes=60,
                        marker="D",
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=2,
                        label=label,
                        hover_data=eqtl_hover_builder.build_dataframe(eqtl_data),
                    )
                    self._backend.add_simple_legend(ax, label, loc="upper right")

            self._backend.set_ylabel(ax, r"$-\log_{10}$ P (eQTL)")
            eqtl_sig_line = -np.log10(eqtl_threshold)
            self._backend.axhline(
                ax,
                y=eqtl_sig_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.65,
            )
            self._backend.hide_spines(ax, ["top", "right"])
            panel_idx += 1

        if genes_df is not None:
            gene_ax = axes[panel_idx]
            plot_gene_track_generic(
                gene_ax, self._backend, genes_df, chrom, start, end, exons_df
            )
            self._backend.hide_spines(gene_ax, ["top", "right", "left"])
            panel_idx += 1

        if heatmap_data is not None:
            heatmap_ax = axes[panel_idx]
            filtered_matrix, x_positions, filtered_snp_ids = heatmap_data
            lead_snp_id = None
            if lead_positions and lead_positions[0] is not None:
                first_gwas = gwas_dfs[0]
                if rs_col in first_gwas.columns:
                    lead_row = first_gwas[first_gwas[pos_col] == lead_positions[0]]
                    if not lead_row.empty:
                        lead_snp_id = lead_row[rs_col].iloc[0]
            self._render_heatmap_panel(
                ax=heatmap_ax,
                fig=fig,
                ld_matrix=filtered_matrix,
                x_positions=x_positions,
                snp_ids=filtered_snp_ids,
                metric=ld_heatmap_metric,
                lead_snp_id=lead_snp_id,
                start=start,
                end=end,
            )
            self._backend.set_xlabel(heatmap_ax, f"Chromosome {chrom} (Mb)")
        elif genes_df is not None:
            self._backend.set_xlabel(gene_ax, f"Chromosome {chrom} (Mb)")
        else:
            self._backend.set_xlabel(axes[-1], f"Chromosome {chrom} (Mb)")

        for ax in axes:
            self._backend.format_xaxis_mb(ax)

        self._backend.finalize_layout(fig, hspace=0.1)

        for ax, texts in all_snp_label_texts:
            self._backend.adjust_snp_labels(ax, texts)

        return fig

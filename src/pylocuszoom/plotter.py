"""Main LocusZoomPlotter class for regional association plots.

Orchestrates all components (LD coloring, gene track, recombination overlay,
SNP labels) into a unified plotting interface.

Supports multiple backends:
- matplotlib (default): Static publication-quality plots
- plotly: Interactive HTML with hover tooltips
- bokeh: Interactive HTML for dashboards
"""

from pathlib import Path
from typing import Any, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

from .backends import BackendType, get_backend
from .backends.hover import HoverConfig, HoverDataBuilder
from .colors import (
    EQTL_NEGATIVE_BINS,
    EQTL_POSITIVE_BINS,
    LD_BINS,
    LEAD_SNP_COLOR,
    PIP_LINE_COLOR,
    get_credible_set_color,
    get_eqtl_color,
    get_ld_bin,
    get_ld_color_palette,
    get_phewas_category_palette,
)
from .config import PlotConfig, StackedPlotConfig
from .ensembl import get_genes_for_region
from .eqtl import validate_eqtl_df
from .finemapping import (
    get_credible_sets,
    prepare_finemapping_for_plotting,
)
from .forest import validate_forest_df
from .gene_track import (
    assign_gene_positions,
    plot_gene_track_generic,
)
from .ld import calculate_ld, find_plink
from .logging import enable_logging, logger
from .manhattan import (
    prepare_categorical_data,
    prepare_manhattan_data,
)
from .phewas import validate_phewas_df
from .qq import prepare_qq_data
from .recombination import (
    RECOMB_COLOR,
    download_canine_recombination_maps,
    get_default_data_dir,
    get_recombination_rate_for_region,
)
from .utils import normalize_chrom, validate_genes_df, validate_gwas_df

# Default significance threshold: 5e-8 (genome-wide significance)
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8
DEFAULT_GENOMEWIDE_LINE = -np.log10(DEFAULT_GENOMEWIDE_THRESHOLD)

# Manhattan/QQ plot styling constants
MANHATTAN_POINT_SIZE = 10
MANHATTAN_CATEGORICAL_POINT_SIZE = 30
QQ_POINT_SIZE = 10  # Match Manhattan point size for consistency
POINT_EDGE_COLOR = "black"
MANHATTAN_EDGE_WIDTH = 0.1
QQ_EDGE_WIDTH = 0.02
SIGNIFICANCE_LINE_COLOR = "red"
QQ_POINT_COLOR = "#1f77b4"
QQ_CI_COLOR = "#CCCCCC"
QQ_CI_ALPHA = 0.5


class LocusZoomPlotter:
    """Regional association plot generator with LD coloring and annotations.

    Creates LocusZoom-style regional plots with:
    - LD coloring based on R² with lead variant
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
        """Initialize the plotter.

        Args:
            species: Species name ('canine', 'feline', or None for custom).
            genome_build: Genome build for coordinate system.
            backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').
            plink_path: Path to PLINK executable for LD calculation.
            recomb_data_dir: Directory containing recombination maps.
            genomewide_threshold: P-value threshold for significance line.
            log_level: Logging level.
            auto_genes: If True, automatically fetch genes from Ensembl when
                genes_df is not provided. Default False for backward compatibility.
        """
        # Configure logging
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
        self._auto_genes = auto_genes

        # Cache for loaded data
        self._recomb_cache = {}

    @staticmethod
    def _default_build(species: str) -> Optional[str]:
        """Get default genome build for species."""
        builds = {"canine": "canfam3.1", "feline": "felCat9"}
        return builds.get(species)

    def _ensure_recomb_maps(self) -> Optional[Path]:
        """Ensure recombination maps are downloaded.

        Returns path to recombination map directory, or None if not available.
        """
        if self.species == "canine":
            if self.recomb_data_dir:
                return Path(self.recomb_data_dir)
            # Check if already downloaded
            default_dir = get_default_data_dir()
            if (
                default_dir.exists()
                and len(list(default_dir.glob("chr*_recomb.tsv"))) >= 39
            ):  # 38 autosomes + X
                return default_dir
            # Download
            try:
                return download_canine_recombination_maps()
            except (requests.RequestException, OSError, IOError) as e:
                # Expected network/file errors - graceful fallback
                logger.warning(f"Could not download recombination maps: {e}")
                return None
            except Exception as e:
                # JUSTIFICATION: Download failure should not prevent plotting.
                # We catch broadly here because graceful degradation is acceptable
                # for optional recombination map downloads. Error-level logging
                # ensures the issue is visible.
                logger.error(f"Unexpected error downloading recombination maps: {e}")
                return None
        elif self.recomb_data_dir:
            return Path(self.recomb_data_dir)
        return None

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
        except FileNotFoundError:
            return None

    def _transform_pvalues(self, df: pd.DataFrame, p_col: str) -> pd.DataFrame:
        """Add neglog10p column with -log10 transformed p-values.

        Clips extremely small p-values to 1e-300 to avoid -inf.

        Args:
            df: DataFrame with p-value column.
            p_col: Name of p-value column.

        Returns:
            DataFrame with neglog10p column added.
        """
        df["neglog10p"] = -np.log10(df[p_col].clip(lower=1e-300))
        return df

    def plot(
        self,
        gwas_df: pd.DataFrame,
        *,
        chrom: int,
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
    ) -> Any:
        """Create a regional association plot.

        Args:
            gwas_df: GWAS results DataFrame.
            chrom: Chromosome number.
            start: Start position in base pairs.
            end: End position in base pairs.
            pos_col: Column name for genomic position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP identifier.
            snp_labels: Whether to show SNP labels on plot.
            label_top_n: Number of top SNPs to label.
            show_recombination: Whether to show recombination rate overlay.
            figsize: Figure size as (width, height) in inches.
            lead_pos: Position of lead SNP to highlight. For stacked plots with
                multiple regions, use plot_stacked() with lead_positions (plural).
            ld_reference_file: Path to PLINK binary fileset for LD calculation.
            ld_col: Column name for pre-computed LD (R^2) values.
            genes_df: Gene annotations with chr, start, end, gene_name.
            exons_df: Exon annotations with chr, start, end, gene_name.
            recomb_df: Pre-loaded recombination rate data.
                If None and show_recombination=True, loads from species default.

        Returns:
            Figure object (type depends on backend).

        Raises:
            ValidationError: If parameters or DataFrame columns are invalid.

        Example:
            >>> fig = plotter.plot(
            ...     gwas_df,
            ...     chrom=1, start=1000000, end=2000000,
            ...     lead_pos=1500000, snp_labels=True,
            ... )
        """
        # Validate parameters via Pydantic
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

        # Validate inputs
        validate_gwas_df(gwas_df, pos_col=pos_col, p_col=p_col)

        # Auto-fetch genes if enabled and not provided
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

        # Prevent auto-display in interactive environments
        plt.ioff()

        # Prepare data
        df = gwas_df.copy()

        # Validate p-values and warn about issues
        p_values = df[p_col]
        nan_count = p_values.isna().sum()
        if nan_count > 0:
            logger.warning(
                f"GWAS data contains {nan_count} NaN p-values which will be excluded"
            )
        invalid_count = ((p_values < 0) | (p_values > 1)).sum()
        if invalid_count > 0:
            logger.warning(
                f"GWAS data contains {invalid_count} p-values outside [0, 1] range"
            )
        clipped_count = (p_values < 1e-300).sum()
        if clipped_count > 0:
            logger.debug(f"Clipping {clipped_count} p-values below 1e-300 to 1e-300")

        df = self._transform_pvalues(df, p_col)

        # Calculate LD if reference file provided
        if ld_reference_file and lead_pos and ld_col is None:
            # Check if rs_col exists before attempting LD calculation
            if rs_col not in df.columns:
                logger.warning(
                    f"Cannot calculate LD: column '{rs_col}' not found in GWAS data. "
                    f"Provide rs_col parameter or add SNP IDs to DataFrame."
                )
            else:
                lead_snp_row = df[df[pos_col] == lead_pos]
                if not lead_snp_row.empty:
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
                        df = df.merge(ld_df, left_on=rs_col, right_on="SNP", how="left")
                        ld_col = "R2"

        # Load recombination data if needed
        if show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(chrom, start, end)

        # Create figure layout
        fig, ax, gene_ax = self._create_figure(genes_df, chrom, start, end, figsize)

        # Plot association data
        self._plot_association(ax, df, pos_col, ld_col, lead_pos, rs_col, p_col)

        # Add significance line
        self._backend.axhline(
            ax,
            y=self._genomewide_line,
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=0.65,
            zorder=1,
        )

        # Add SNP labels (capability check - interactive backends use hover tooltips)
        if snp_labels and rs_col in df.columns and label_top_n > 0 and not df.empty:
            if self._backend.supports_snp_labels:
                self._backend.add_snp_labels(
                    ax,
                    df,
                    pos_col=pos_col,
                    neglog10p_col="neglog10p",
                    rs_col=rs_col,
                    label_top_n=label_top_n,
                    genes_df=genes_df,
                    chrom=chrom,
                )

        # Add recombination overlay (all backends with secondary axis support)
        if recomb_df is not None and not recomb_df.empty:
            if self._backend.supports_secondary_axis:
                self._add_recombination_overlay(ax, recomb_df, start, end)

        # Format axes
        self._backend.set_ylabel(ax, r"$-\log_{10}$ P")
        self._backend.set_xlim(ax, start, end)
        # When recombination overlay is present, keep right spine for secondary y-axis
        has_recomb = recomb_df is not None and not recomb_df.empty
        if has_recomb and self._backend.supports_secondary_axis:
            self._backend.hide_spines(ax, ["top"])
        else:
            self._backend.hide_spines(ax, ["top", "right"])

        # Add LD legend (all backends)
        if ld_col is not None and ld_col in df.columns:
            self._backend.add_ld_legend(ax, LD_BINS, LEAD_SNP_COLOR)

        # Plot gene track (all backends use generic function)
        if genes_df is not None and gene_ax is not None:
            plot_gene_track_generic(
                gene_ax, self._backend, genes_df, chrom, start, end, exons_df
            )
            self._backend.set_xlabel(gene_ax, f"Chromosome {chrom} (Mb)")
            self._backend.hide_spines(gene_ax, ["top", "right", "left"])
            # Format both axes for interactive backends (they don't share x-axis)
            self._backend.format_xaxis_mb(gene_ax)
        else:
            self._backend.set_xlabel(ax, f"Chromosome {chrom} (Mb)")

        # Format x-axis with Mb labels (association axis always needs formatting)
        self._backend.format_xaxis_mb(ax)

        # Adjust layout
        self._backend.finalize_layout(fig, hspace=0.1)

        return fig

    def _create_figure(
        self,
        genes_df: Optional[pd.DataFrame],
        chrom: int,
        start: int,
        end: int,
        figsize: Tuple[int, int],
    ) -> Tuple[Any, Any, Optional[Any]]:
        """Create figure with optional gene track."""
        if genes_df is not None:
            # Calculate dynamic height based on gene rows
            chrom_str = normalize_chrom(chrom)
            region_genes = genes_df[
                (
                    genes_df["chr"].astype(str).str.replace("chr", "", regex=False)
                    == chrom_str
                )
                & (genes_df["end"] >= start)
                & (genes_df["start"] <= end)
            ]
            if not region_genes.empty:
                temp_positions = assign_gene_positions(
                    region_genes.sort_values("start"), start, end
                )
                n_gene_rows = max(temp_positions) + 1 if temp_positions else 1
            else:
                n_gene_rows = 1

            base_gene_height = 1.0
            per_row_height = 0.5
            gene_track_height = base_gene_height + (n_gene_rows - 1) * per_row_height
            assoc_height = figsize[1] * 0.6
            total_height = assoc_height + gene_track_height

            fig, axes = self._backend.create_figure(
                n_panels=2,
                height_ratios=[assoc_height, gene_track_height],
                figsize=(figsize[0], total_height),
                sharex=True,
            )
            return fig, axes[0], axes[1]
        else:
            fig, axes = self._backend.create_figure(
                n_panels=1,
                height_ratios=[1.0],
                figsize=(figsize[0], figsize[1] * 0.75),
            )
            return fig, axes[0], None

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
        # Build hover data using HoverDataBuilder
        hover_config = HoverConfig(
            snp_col=rs_col if rs_col and rs_col in df.columns else None,
            pos_col=pos_col if pos_col in df.columns else None,
            p_col=p_col if p_col and p_col in df.columns else None,
            ld_col=ld_col if ld_col and ld_col in df.columns else None,
        )
        hover_builder = HoverDataBuilder(hover_config)

        # LD-based coloring
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
            # Default: grey points
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

        # Highlight lead SNP with larger, more prominent marker
        if lead_pos is not None:
            lead_snp = df[df[pos_col] == lead_pos]
            if not lead_snp.empty:
                self._backend.scatter(
                    ax,
                    lead_snp[pos_col],
                    lead_snp["neglog10p"],
                    colors=LEAD_SNP_COLOR,
                    sizes=120,  # Larger than regular points for visibility
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
        """Add recombination overlay for all backends.

        Creates a secondary y-axis with recombination rate line and fill.
        Uses backend-agnostic secondary axis methods that work across
        matplotlib, plotly, and bokeh.
        """
        # Filter to region
        region_recomb = recomb_df[
            (recomb_df["pos"] >= start) & (recomb_df["pos"] <= end)
        ].copy()

        if region_recomb.empty:
            return

        # Create secondary y-axis
        twin_result = self._backend.create_twin_axis(ax)

        # Matplotlib returns the twin Axes object itself - use it for drawing
        # Plotly returns tuple (fig, row, secondary_y_name)
        # Bokeh returns string "secondary"
        from matplotlib.axes import Axes

        if isinstance(twin_result, Axes):
            # Matplotlib: use the twin axis for all secondary axis operations
            secondary_ax = twin_result
            secondary_y = None  # Not used for matplotlib
        elif isinstance(twin_result, tuple):
            # Plotly: use original ax, specify y-axis via yaxis_name
            secondary_ax = ax
            _, _, secondary_y = twin_result
        else:
            # Bokeh: use original ax, specify y-axis via yaxis_name
            secondary_ax = ax
            secondary_y = twin_result

        # Plot fill under curve
        self._backend.fill_between_secondary(
            secondary_ax,
            region_recomb["pos"],
            0,
            region_recomb["rate"],
            color=RECOMB_COLOR,
            alpha=0.15,
            yaxis_name=secondary_y,
        )

        # Plot recombination rate line
        self._backend.line_secondary(
            secondary_ax,
            region_recomb["pos"],
            region_recomb["rate"],
            color=RECOMB_COLOR,
            linewidth=2.5,
            alpha=0.8,
            yaxis_name=secondary_y,
        )

        # Set y-axis limits and label - scale to fit data with headroom
        max_rate = region_recomb["rate"].max()
        self._backend.set_secondary_ylim(
            secondary_ax, 0, max(max_rate * 1.3, 10), yaxis_name=secondary_y
        )
        self._backend.set_secondary_ylabel(
            secondary_ax,
            "Recombination rate (cM/Mb)",
            color="black",  # Use black for readability (line/fill color remains light blue)
            fontsize=9,
            yaxis_name=secondary_y,
        )

        # Hide top spine on the secondary axis (matplotlib twin axis has its own frame)
        if isinstance(twin_result, Axes):
            secondary_ax.spines["top"].set_visible(False)

    def _plot_finemapping(
        self,
        ax: Any,
        df: pd.DataFrame,
        pos_col: str = "pos",
        pip_col: str = "pip",
        cs_col: Optional[str] = "cs",
        show_credible_sets: bool = True,
        pip_threshold: float = 0.0,
    ) -> None:
        """Plot fine-mapping results (PIP line with credible set coloring).

        Args:
            ax: Matplotlib axes object.
            df: Fine-mapping DataFrame with pos and pip columns.
            pos_col: Column name for position.
            pip_col: Column name for posterior inclusion probability.
            cs_col: Column name for credible set assignment (optional).
            show_credible_sets: Whether to color points by credible set.
            pip_threshold: Minimum PIP to display as scatter point.
        """
        # Build hover data using HoverDataBuilder
        extra_cols = {pip_col: "PIP"}
        if cs_col and cs_col in df.columns:
            extra_cols[cs_col] = "Credible Set"
        hover_config = HoverConfig(
            pos_col=pos_col if pos_col in df.columns else None,
            extra_cols=extra_cols,
        )
        hover_builder = HoverDataBuilder(hover_config)

        # Sort by position for line plotting
        df = df.sort_values(pos_col)

        # Plot PIP as line
        self._backend.line(
            ax,
            df[pos_col],
            df[pip_col],
            color=PIP_LINE_COLOR,
            linewidth=1.5,
            alpha=0.8,
            zorder=1,
        )

        # Check if credible sets are available
        has_cs = cs_col is not None and cs_col in df.columns and show_credible_sets
        credible_sets = get_credible_sets(df, cs_col) if has_cs else []

        if credible_sets:
            # Plot points colored by credible set
            for cs_id in credible_sets:
                cs_data = df[df[cs_col] == cs_id]
                color = get_credible_set_color(cs_id)
                self._backend.scatter(
                    ax,
                    cs_data[pos_col],
                    cs_data[pip_col],
                    colors=color,
                    sizes=50,
                    marker="o",
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=3,
                    hover_data=hover_builder.build_dataframe(cs_data),
                )
            # Plot variants not in any credible set
            non_cs_data = df[(df[cs_col].isna()) | (df[cs_col] == 0)]
            if not non_cs_data.empty and pip_threshold > 0:
                non_cs_data = non_cs_data[non_cs_data[pip_col] >= pip_threshold]
                if not non_cs_data.empty:
                    self._backend.scatter(
                        ax,
                        non_cs_data[pos_col],
                        non_cs_data[pip_col],
                        colors="#BEBEBE",
                        sizes=30,
                        marker="o",
                        edgecolor="black",
                        linewidth=0.3,
                        zorder=2,
                        hover_data=hover_builder.build_dataframe(non_cs_data),
                    )
        else:
            # No credible sets - show all points above threshold
            if pip_threshold > 0:
                high_pip = df[df[pip_col] >= pip_threshold]
                if not high_pip.empty:
                    self._backend.scatter(
                        ax,
                        high_pip[pos_col],
                        high_pip[pip_col],
                        colors=PIP_LINE_COLOR,
                        sizes=50,
                        marker="o",
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=3,
                        hover_data=hover_builder.build_dataframe(high_pip),
                    )

    def plot_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        *,
        chrom: int,
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
        finemapping_df: Optional[pd.DataFrame] = None,
        finemapping_cs_col: Optional[str] = "cs",
        recomb_df: Optional[pd.DataFrame] = None,
    ) -> Any:
        """Create stacked regional association plots for multiple GWAS.

        Vertically stacks multiple GWAS results for comparison, with shared
        x-axis and optional gene track at the bottom.

        Args:
            gwas_dfs: List of GWAS results DataFrames to stack.
            chrom: Chromosome number.
            start: Start position in base pairs.
            end: End position in base pairs.
            pos_col: Column name for genomic position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP identifier.
            snp_labels: Whether to show SNP labels on plot.
            label_top_n: Number of top SNPs to label (default 3 for stacked).
            show_recombination: Whether to show recombination rate overlay.
            figsize: Figure size as (width, height) in inches.
            ld_reference_file: Single PLINK fileset (broadcast to all panels).
            ld_col: Column name for pre-computed LD (R^2) values.
            lead_positions: List of lead SNP positions, one per region. For single
                region plots, use plot() with lead_pos (singular).
            panel_labels: List of panel labels (one per panel).
            ld_reference_files: List of PLINK filesets (one per panel).
            genes_df: Gene annotations for bottom track.
            exons_df: Exon annotations for gene track.
            eqtl_df: eQTL data to display as additional panel.
            eqtl_gene: Filter eQTL data to this target gene.
            finemapping_df: Fine-mapping/SuSiE results with pos and pip columns.
                Displayed as PIP line with optional credible set coloring.
            finemapping_cs_col: Column name for credible set assignment.
            recomb_df: Pre-loaded recombination rate data.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_stacked(
            ...     [gwas_height, gwas_bmi, gwas_whr],
            ...     chrom=1, start=1000000, end=2000000,
            ...     panel_labels=["Height", "BMI", "WHR"],
            ... )
        """
        # Validate parameters via Pydantic
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

        # Validate list lengths match
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

        # Validate inputs
        for i, df in enumerate(gwas_dfs):
            validate_gwas_df(df, pos_col=pos_col, p_col=p_col)
        if genes_df is not None:
            validate_genes_df(genes_df)
        if eqtl_df is not None:
            validate_eqtl_df(eqtl_df)

        # Handle lead positions
        if lead_positions is None:
            lead_positions = []
            for df in gwas_dfs:
                region_df = df[(df[pos_col] >= start) & (df[pos_col] <= end)]
                if not region_df.empty:
                    # Filter out NaN p-values for lead SNP detection
                    valid_p = region_df[p_col].dropna()
                    if valid_p.empty:
                        logger.warning(
                            "All p-values in region are NaN, cannot determine lead SNP"
                        )
                        lead_positions.append(None)
                    else:
                        lead_idx = valid_p.idxmin()
                        lead_positions.append(int(region_df.loc[lead_idx, pos_col]))
                else:
                    lead_positions.append(None)

        # Handle LD reference files
        if ld_reference_files is None and ld_reference_file is not None:
            ld_reference_files = [ld_reference_file] * n_gwas

        # Calculate panel layout
        panel_height = 2.5  # inches per GWAS panel
        eqtl_height = 2.0 if eqtl_df is not None else 0
        finemapping_height = 1.5 if finemapping_df is not None else 0

        # Gene track height
        if genes_df is not None:
            chrom_str = normalize_chrom(chrom)
            region_genes = genes_df[
                (
                    genes_df["chr"].astype(str).str.replace("chr", "", regex=False)
                    == chrom_str
                )
                & (genes_df["end"] >= start)
                & (genes_df["start"] <= end)
            ]
            if not region_genes.empty:
                temp_positions = assign_gene_positions(
                    region_genes.sort_values("start"), start, end
                )
                n_gene_rows = max(temp_positions) + 1 if temp_positions else 1
            else:
                n_gene_rows = 1
            gene_track_height = 1.0 + (n_gene_rows - 1) * 0.5
        else:
            gene_track_height = 0

        # Calculate total panels and heights
        n_panels = (
            n_gwas
            + (1 if finemapping_df is not None else 0)
            + (1 if eqtl_df is not None else 0)
            + (1 if genes_df is not None else 0)
        )
        height_ratios = [panel_height] * n_gwas
        if finemapping_df is not None:
            height_ratios.append(finemapping_height)
        if eqtl_df is not None:
            height_ratios.append(eqtl_height)
        if genes_df is not None:
            height_ratios.append(gene_track_height)

        # Calculate figure height
        total_height = figsize[1] if figsize[1] else sum(height_ratios)
        actual_figsize = (figsize[0], total_height)

        logger.debug(
            f"Creating stacked plot with {n_panels} panels for chr{chrom}:{start}-{end}"
        )

        # Load recombination data if needed
        if show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(chrom, start, end)

        # Create figure using backend
        fig, axes = self._backend.create_figure(
            n_panels=n_panels,
            height_ratios=height_ratios,
            figsize=actual_figsize,
            sharex=True,
        )

        # Plot each GWAS panel
        for i, (gwas_df, lead_pos) in enumerate(zip(gwas_dfs, lead_positions)):
            ax = axes[i]
            df = gwas_df.copy()
            df = self._transform_pvalues(df, p_col)

            # Use pre-computed LD or calculate from reference
            panel_ld_col = ld_col
            if ld_reference_files and ld_reference_files[i] and lead_pos and not ld_col:
                lead_snp_row = df[df[pos_col] == lead_pos]
                if not lead_snp_row.empty and rs_col in df.columns:
                    lead_snp_id = lead_snp_row[rs_col].iloc[0]
                    ld_df = calculate_ld(
                        bfile_path=ld_reference_files[i],
                        lead_snp=lead_snp_id,
                        window_kb=max((end - start) // 1000, 500),
                        plink_path=self.plink_path,
                        species=self.species,
                    )
                    if not ld_df.empty:
                        df = df.merge(ld_df, left_on=rs_col, right_on="SNP", how="left")
                        panel_ld_col = "R2"

            # Plot association
            self._plot_association(
                ax, df, pos_col, panel_ld_col, lead_pos, rs_col, p_col
            )

            # Add significance line
            self._backend.axhline(
                ax,
                y=self._genomewide_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.65,
                zorder=1,
            )

            # Add SNP labels (capability check - interactive backends use hover tooltips)
            if snp_labels and rs_col in df.columns and label_top_n > 0 and not df.empty:
                if self._backend.supports_snp_labels:
                    self._backend.add_snp_labels(
                        ax,
                        df,
                        pos_col=pos_col,
                        neglog10p_col="neglog10p",
                        rs_col=rs_col,
                        label_top_n=label_top_n,
                        genes_df=genes_df,
                        chrom=chrom,
                    )

            # Add recombination overlay (only on first panel, all backends)
            if i == 0 and recomb_df is not None and not recomb_df.empty:
                if self._backend.supports_secondary_axis:
                    self._add_recombination_overlay(ax, recomb_df, start, end)

            # Format axes
            self._backend.set_ylabel(ax, r"$-\log_{10}$ P")
            self._backend.set_xlim(ax, start, end)
            self._backend.hide_spines(ax, ["top", "right"])

            # Add panel label
            if panel_labels and i < len(panel_labels):
                self._backend.add_panel_label(ax, panel_labels[i])

            # Add LD legend (only on first panel, all backends)
            if i == 0 and panel_ld_col is not None and panel_ld_col in df.columns:
                self._backend.add_ld_legend(ax, LD_BINS, LEAD_SNP_COLOR)

        # Track current panel index
        panel_idx = n_gwas

        # Plot fine-mapping panel if provided
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
                self._plot_finemapping(
                    ax,
                    fm_data,
                    pos_col="pos",
                    pip_col="pip",
                    cs_col=finemapping_cs_col,
                    show_credible_sets=True,
                    pip_threshold=0.01,
                )

                # Add legend for credible sets (all backends)
                credible_sets = get_credible_sets(fm_data, finemapping_cs_col)
                if credible_sets:
                    self._backend.add_finemapping_legend(
                        ax, credible_sets, get_credible_set_color
                    )

            self._backend.set_ylabel(ax, "PIP")
            self._backend.set_ylim(ax, -0.05, 1.05)
            self._backend.hide_spines(ax, ["top", "right"])
            panel_idx += 1

        # Plot eQTL panel if provided
        eqtl_panel_idx = panel_idx
        if eqtl_df is not None:
            ax = axes[eqtl_panel_idx]
            eqtl_data = eqtl_df.copy()

            # Filter by gene if specified
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

            # Filter by region (position and chromosome)
            if "pos" in eqtl_data.columns:
                mask = (eqtl_data["pos"] >= start) & (eqtl_data["pos"] <= end)
                # Also filter by chromosome if column exists
                if "chr" in eqtl_data.columns:
                    chrom_str = str(chrom).replace("chr", "")
                    eqtl_chrom = (
                        eqtl_data["chr"].astype(str).str.replace("chr", "", regex=False)
                    )
                    mask = mask & (eqtl_chrom == chrom_str)
                eqtl_data = eqtl_data[mask]

            if not eqtl_data.empty:
                eqtl_data = self._transform_pvalues(eqtl_data, "p_value")

                # Build hover data using HoverDataBuilder
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

                # Check if effect_size column exists for directional coloring
                has_effect = "effect_size" in eqtl_data.columns

                if has_effect:
                    # Vectorized plotting: split by sign, assign colors in bulk
                    pos_effects = eqtl_data[eqtl_data["effect_size"] >= 0]
                    neg_effects = eqtl_data[eqtl_data["effect_size"] < 0]

                    # Vectorized color assignment using apply
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

                    # Add eQTL effect legend (all backends)
                    self._backend.add_eqtl_legend(
                        ax, EQTL_POSITIVE_BINS, EQTL_NEGATIVE_BINS
                    )
                else:
                    # No effect sizes - plot as diamonds
                    # Only show gene in label if filtering was actually applied
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
            self._backend.axhline(
                ax,
                y=self._genomewide_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.65,
            )
            self._backend.hide_spines(ax, ["top", "right"])
            panel_idx += 1

        # Plot gene track (all backends use generic function)
        if genes_df is not None:
            gene_ax = axes[panel_idx]
            plot_gene_track_generic(
                gene_ax, self._backend, genes_df, chrom, start, end, exons_df
            )
            self._backend.set_xlabel(gene_ax, f"Chromosome {chrom} (Mb)")
            self._backend.hide_spines(gene_ax, ["top", "right", "left"])
        else:
            # Set x-label on bottom panel
            self._backend.set_xlabel(axes[-1], f"Chromosome {chrom} (Mb)")

        # Format x-axis (call for all axes - Plotly needs each subplot formatted)
        for ax in axes:
            self._backend.format_xaxis_mb(ax)

        # Adjust layout
        self._backend.finalize_layout(fig, hspace=0.1)

        return fig

    def plot_phewas(
        self,
        phewas_df: pd.DataFrame,
        variant_id: str,
        phenotype_col: str = "phenotype",
        p_col: str = "p_value",
        category_col: str = "category",
        effect_col: Optional[str] = None,
        significance_threshold: float = 5e-8,
        figsize: Tuple[float, float] = (10, 8),
    ) -> Any:
        """Create a PheWAS (Phenome-Wide Association Study) plot.

        Shows associations of a single variant across multiple phenotypes,
        with phenotypes grouped by category and colored accordingly.

        Args:
            phewas_df: DataFrame with phenotype associations.
            variant_id: Variant identifier (e.g., "rs12345") for plot title.
            phenotype_col: Column name for phenotype names.
            p_col: Column name for p-values.
            category_col: Column name for phenotype categories.
            effect_col: Optional column name for effect direction (beta/OR).
            significance_threshold: P-value threshold for significance line.
            figsize: Figure size as (width, height).

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_phewas(
            ...     phewas_df,
            ...     variant_id="rs12345",
            ...     category_col="category",
            ... )
        """
        validate_phewas_df(phewas_df, phenotype_col, p_col, category_col)

        df = phewas_df.copy()
        df = self._transform_pvalues(df, p_col)

        # Sort by category then by p-value for consistent ordering
        if category_col in df.columns:
            df = df.sort_values([category_col, p_col])
            categories = df[category_col].unique().tolist()
            palette = get_phewas_category_palette(categories)
        else:
            df = df.sort_values(p_col)
            categories = []
            palette = {}

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Assign y-positions (one per phenotype)
        df["y_pos"] = range(len(df))

        # Plot points by category
        if categories:
            for cat in categories:
                # Handle NaN category: NaN == NaN is False in pandas
                if pd.isna(cat):
                    cat_data = df[df[category_col].isna()]
                else:
                    cat_data = df[df[category_col] == cat]
                # Use upward triangles for positive effects, circles otherwise
                if effect_col and effect_col in cat_data.columns:
                    # Vectorized: split by effect sign, 2 scatter calls per category
                    pos_data = cat_data[cat_data[effect_col] >= 0]
                    neg_data = cat_data[cat_data[effect_col] < 0]

                    if not pos_data.empty:
                        self._backend.scatter(
                            ax,
                            pos_data["neglog10p"],
                            pos_data["y_pos"],
                            colors=palette[cat],
                            sizes=60,
                            marker="^",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
                    if not neg_data.empty:
                        self._backend.scatter(
                            ax,
                            neg_data["neglog10p"],
                            neg_data["y_pos"],
                            colors=palette[cat],
                            sizes=60,
                            marker="v",
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
                else:
                    self._backend.scatter(
                        ax,
                        cat_data["neglog10p"],
                        cat_data["y_pos"],
                        colors=palette[cat],
                        sizes=60,
                        marker="o",
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=2,
                    )
        else:
            self._backend.scatter(
                ax,
                df["neglog10p"],
                df["y_pos"],
                colors="#4169E1",
                sizes=60,
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
            )

        # Add significance threshold line
        sig_line = -np.log10(significance_threshold)
        self._backend.axvline(
            ax, x=sig_line, color="red", linestyle="--", linewidth=1, alpha=0.7
        )

        # Set axis labels and limits
        self._backend.set_xlabel(ax, r"$-\log_{10}$ P")
        self._backend.set_ylabel(ax, "Phenotype")
        self._backend.set_ylim(ax, -0.5, len(df) - 0.5)

        # Set y-tick labels to phenotype names
        self._backend.set_yticks(
            ax,
            positions=df["y_pos"].tolist(),
            labels=df[phenotype_col].tolist(),
            fontsize=8,
        )

        self._backend.set_title(ax, f"PheWAS: {variant_id}")
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

    def plot_forest(
        self,
        forest_df: pd.DataFrame,
        variant_id: str,
        study_col: str = "study",
        effect_col: str = "effect",
        ci_lower_col: str = "ci_lower",
        ci_upper_col: str = "ci_upper",
        weight_col: Optional[str] = None,
        null_value: float = 0.0,
        effect_label: str = "Effect Size",
        figsize: Tuple[float, float] = (8, 6),
    ) -> Any:
        """Create a forest plot showing effect sizes with confidence intervals.

        Args:
            forest_df: DataFrame with effect sizes and confidence intervals.
            variant_id: Variant identifier for plot title.
            study_col: Column name for study/phenotype names.
            effect_col: Column name for effect sizes.
            ci_lower_col: Column name for lower confidence interval.
            ci_upper_col: Column name for upper confidence interval.
            weight_col: Optional column for study weights (affects marker size).
            null_value: Reference value for null effect (0 for beta, 1 for OR).
            effect_label: X-axis label.
            figsize: Figure size as (width, height).

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_forest(
            ...     forest_df,
            ...     variant_id="rs12345",
            ...     effect_label="Odds Ratio",
            ...     null_value=1.0,
            ... )
        """
        validate_forest_df(forest_df, study_col, effect_col, ci_lower_col, ci_upper_col)

        df = forest_df.copy()

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Assign y-positions (reverse so first study is at top)
        df["y_pos"] = range(len(df) - 1, -1, -1)

        # Calculate marker sizes from weights
        if weight_col and weight_col in df.columns:
            # Scale weights to marker sizes (min 40, max 200)
            weights = df[weight_col]
            min_size, max_size = 40, 200
            weight_range = weights.max() - weights.min()
            if weight_range > 0:
                sizes = min_size + (weights - weights.min()) / weight_range * (
                    max_size - min_size
                )
            else:
                sizes = (min_size + max_size) / 2
        else:
            sizes = 80

        # Calculate error bar extents
        xerr_lower = df[effect_col] - df[ci_lower_col]
        xerr_upper = df[ci_upper_col] - df[effect_col]

        # Plot error bars (confidence intervals)
        self._backend.errorbar_h(
            ax,
            x=df[effect_col],
            y=df["y_pos"],
            xerr_lower=xerr_lower,
            xerr_upper=xerr_upper,
            color="black",
            linewidth=1.5,
            capsize=3,
            zorder=2,
        )

        # Plot effect size markers
        self._backend.scatter(
            ax,
            df[effect_col],
            df["y_pos"],
            colors="#4169E1",
            sizes=sizes,
            marker="s",  # square markers typical for forest plots
            edgecolor="black",
            linewidth=0.5,
            zorder=3,
        )

        # Add null effect line
        self._backend.axvline(
            ax, x=null_value, color="grey", linestyle="--", linewidth=1, alpha=0.7
        )

        # Set axis labels and limits
        self._backend.set_xlabel(ax, effect_label)
        self._backend.set_ylim(ax, -0.5, len(df) - 0.5)

        # Ensure x-axis includes the null value with some padding
        x_min = min(df[ci_lower_col].min(), null_value)
        x_max = max(df[ci_upper_col].max(), null_value)
        x_padding = (x_max - x_min) * 0.1
        self._backend.set_xlim(ax, x_min - x_padding, x_max + x_padding)

        # Set y-tick labels to study names
        self._backend.set_yticks(
            ax,
            positions=df["y_pos"].tolist(),
            labels=df[study_col].tolist(),
            fontsize=10,
        )

        self._backend.set_title(ax, f"Forest Plot: {variant_id}")
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

    def plot_manhattan(
        self,
        df: pd.DataFrame,
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        category_col: Optional[str] = None,
        category_order: Optional[List[str]] = None,
        significance_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        figsize: Tuple[float, float] = (12, 5),
        title: Optional[str] = None,
    ) -> Any:
        """Create a Manhattan plot.

        Shows associations across the genome with points colored by chromosome.
        Supports both standard Manhattan plots (genomic positions) and
        categorical Manhattan plots (PheWAS-style).

        Args:
            df: DataFrame with GWAS results.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            custom_chrom_order: Custom chromosome order (overrides species).
            category_col: If provided, creates a categorical Manhattan plot
                (like PheWAS) using this column instead of genomic positions.
            category_order: Custom category order for categorical plots.
            significance_threshold: P-value threshold for genome-wide significance
                line. Set to None to disable.
            figsize: Figure size as (width, height).
            title: Plot title. Defaults to "Manhattan Plot".

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> # Standard Manhattan plot
            >>> fig = plotter.plot_manhattan(gwas_df, species="human")
            >>>
            >>> # Categorical Manhattan (PheWAS-style)
            >>> fig = plotter.plot_manhattan(
            ...     phewas_df,
            ...     category_col="phenotype_category",
            ...     p_col="pvalue",
            ... )
        """
        # Categorical Manhattan plot
        if category_col is not None:
            return self._plot_manhattan_categorical(
                df=df,
                category_col=category_col,
                p_col=p_col,
                category_order=category_order,
                significance_threshold=significance_threshold,
                figsize=figsize,
                title=title,
            )

        # Standard Manhattan plot
        # Prepare data with cumulative positions and colors
        prepared_df = prepare_manhattan_data(
            df=df,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Plot points and significance line
        chrom_order = prepared_df.attrs["chrom_order"]
        self._render_manhattan_points(ax, prepared_df, chrom_order)
        self._add_significance_line(ax, significance_threshold)

        # Set x-axis ticks to chromosome centers
        chrom_centers = prepared_df.attrs["chrom_centers"]
        positions = [
            chrom_centers[chrom] for chrom in chrom_order if chrom in chrom_centers
        ]
        labels = [chrom for chrom in chrom_order if chrom in chrom_centers]
        self._backend.set_xticks(ax, positions, labels, fontsize=8)

        # Set limits
        x_min = prepared_df["_cumulative_pos"].min()
        x_max = prepared_df["_cumulative_pos"].max()
        x_padding = (x_max - x_min) * 0.01
        self._backend.set_xlim(ax, x_min - x_padding, x_max + x_padding)

        y_max = prepared_df["_neg_log_p"].max()
        self._backend.set_ylim(ax, 0, y_max * 1.1)

        # Labels and title
        self._backend.set_xlabel(ax, "Chromosome", fontsize=12)
        self._backend.set_ylabel(ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_title(ax, title or "Manhattan Plot", fontsize=14)
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

    def _render_manhattan_points(
        self,
        ax: Any,
        prepared_df: pd.DataFrame,
        chrom_order: List[str],
        point_size: int = MANHATTAN_POINT_SIZE,
    ) -> None:
        """Render Manhattan plot scatter points grouped by chromosome.

        Args:
            ax: Axes object from backend.
            prepared_df: DataFrame with _chrom_str, _cumulative_pos, _neg_log_p, _color.
            chrom_order: List of chromosome names in display order.
            point_size: Size of scatter points.
        """
        for chrom in chrom_order:
            chrom_data = prepared_df[prepared_df["_chrom_str"] == chrom]
            if len(chrom_data) > 0:
                self._backend.scatter(
                    ax,
                    chrom_data["_cumulative_pos"],
                    chrom_data["_neg_log_p"],
                    colors=chrom_data["_color"].iloc[0],
                    sizes=point_size,
                    marker="o",
                    edgecolor=POINT_EDGE_COLOR,
                    linewidth=MANHATTAN_EDGE_WIDTH,
                    zorder=2,
                )

    def _add_significance_line(
        self,
        ax: Any,
        threshold: Optional[float],
    ) -> None:
        """Add genome-wide significance threshold line.

        Args:
            ax: Axes object from backend.
            threshold: P-value threshold (e.g., 5e-8). None to skip.
        """
        if threshold is None:
            return
        threshold_line = -np.log10(threshold)
        self._backend.axhline(
            ax,
            y=threshold_line,
            color=SIGNIFICANCE_LINE_COLOR,
            linestyle="--",
            linewidth=1,
            zorder=1,
        )

    def _render_qq_plot(
        self,
        ax: Any,
        qq_df: pd.DataFrame,
        show_confidence_band: bool = True,
    ) -> None:
        """Render QQ plot elements on axes.

        Args:
            ax: Axes object from backend.
            qq_df: Prepared QQ DataFrame with _expected, _observed, _ci_lower, _ci_upper.
            show_confidence_band: Whether to show 95% confidence band.
        """
        if show_confidence_band:
            self._backend.fill_between(
                ax,
                x=qq_df["_expected"],
                y1=qq_df["_ci_lower"],
                y2=qq_df["_ci_upper"],
                color=QQ_CI_COLOR,
                alpha=QQ_CI_ALPHA,
                zorder=1,
            )

        max_val = max(qq_df["_expected"].max(), qq_df["_observed"].max())

        # Diagonal reference line
        self._backend.line(
            ax,
            x=pd.Series([0, max_val]),
            y=pd.Series([0, max_val]),
            color=SIGNIFICANCE_LINE_COLOR,
            linestyle="--",
            linewidth=1,
            zorder=2,
        )

        # QQ points
        self._backend.scatter(
            ax,
            qq_df["_expected"],
            qq_df["_observed"],
            colors=QQ_POINT_COLOR,
            sizes=QQ_POINT_SIZE,
            marker="o",
            edgecolor=POINT_EDGE_COLOR,
            linewidth=QQ_EDGE_WIDTH,
            zorder=3,
        )

        # Set limits
        self._backend.set_xlim(ax, 0, max_val * 1.05)
        self._backend.set_ylim(ax, 0, max_val * 1.05)

    def _plot_manhattan_categorical(
        self,
        df: pd.DataFrame,
        category_col: str,
        p_col: str = "p",
        category_order: Optional[List[str]] = None,
        significance_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        figsize: Tuple[float, float] = (12, 5),
        title: Optional[str] = None,
    ) -> Any:
        """Create a categorical Manhattan plot (PheWAS-style).

        Internal method called by plot_manhattan when category_col is provided.
        """
        # Prepare data
        prepared_df = prepare_categorical_data(
            df=df,
            category_col=category_col,
            p_col=p_col,
            category_order=category_order,
        )

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Plot points by category
        cat_order = prepared_df.attrs["category_order"]
        for cat in cat_order:
            cat_data = prepared_df[prepared_df[category_col] == cat]
            if len(cat_data) > 0:
                self._backend.scatter(
                    ax,
                    cat_data["_x_pos"],
                    cat_data["_neg_log_p"],
                    colors=cat_data["_color"].iloc[0],
                    sizes=MANHATTAN_CATEGORICAL_POINT_SIZE,
                    marker="o",
                    edgecolor=POINT_EDGE_COLOR,
                    linewidth=MANHATTAN_EDGE_WIDTH,
                    zorder=2,
                )

        self._add_significance_line(ax, significance_threshold)

        # Set x-axis ticks
        cat_centers = prepared_df.attrs["category_centers"]
        positions = [cat_centers[cat] for cat in cat_order]
        self._backend.set_xticks(
            ax, positions, cat_order, fontsize=10, rotation=45, ha="right"
        )

        # Set limits
        self._backend.set_xlim(ax, -0.5, len(cat_order) - 0.5)

        y_max = prepared_df["_neg_log_p"].max()
        self._backend.set_ylim(ax, 0, y_max * 1.1)

        # Labels and title
        self._backend.set_xlabel(ax, "Category", fontsize=12)
        self._backend.set_ylabel(ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_title(ax, title or "Categorical Manhattan Plot", fontsize=14)
        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

    def plot_qq(
        self,
        df: pd.DataFrame,
        p_col: str = "p",
        show_confidence_band: bool = True,
        show_lambda: bool = True,
        figsize: Tuple[float, float] = (6, 6),
        title: Optional[str] = None,
    ) -> Any:
        """Create a QQ (quantile-quantile) plot.

        Shows observed vs expected -log10(p) distribution with optional
        95% confidence band and genomic inflation factor (lambda).

        Args:
            df: DataFrame with p-values.
            p_col: Column name for p-value.
            show_confidence_band: If True, show 95% confidence band.
            show_lambda: If True, show genomic inflation factor in title.
            figsize: Figure size as (width, height).
            title: Plot title. If None and show_lambda is True, shows lambda.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_qq(gwas_df, p_col="pvalue")
        """
        # Prepare data
        prepared_df = prepare_qq_data(df, p_col=p_col)

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]

        # Render QQ plot elements
        self._render_qq_plot(ax, prepared_df, show_confidence_band)

        # Labels
        self._backend.set_xlabel(ax, r"Expected $-\log_{10}(p)$", fontsize=12)
        self._backend.set_ylabel(ax, r"Observed $-\log_{10}(p)$", fontsize=12)

        # Title with lambda
        if title:
            plot_title = title
        elif show_lambda:
            lambda_gc = prepared_df.attrs["lambda_gc"]
            plot_title = f"QQ Plot (λ = {lambda_gc:.3f})"
        else:
            plot_title = "QQ Plot"
        self._backend.set_title(ax, plot_title, fontsize=14)

        self._backend.hide_spines(ax, ["top", "right"])
        self._backend.finalize_layout(fig)

        return fig

    def plot_manhattan_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        significance_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        panel_labels: Optional[List[str]] = None,
        figsize: Tuple[float, float] = (12, 8),
        title: Optional[str] = None,
    ) -> Any:
        """Create stacked Manhattan plots for multiple GWAS datasets.

        Vertically stacks multiple Manhattan plots for easy comparison across
        studies or phenotypes.

        Args:
            gwas_dfs: List of GWAS results DataFrames.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            custom_chrom_order: Custom chromosome order (overrides species).
            significance_threshold: P-value threshold for genome-wide significance
                line. Set to None to disable.
            panel_labels: Labels for each panel (one per DataFrame).
            figsize: Figure size as (width, height).
            title: Overall plot title.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_manhattan_stacked(
            ...     [gwas1, gwas2, gwas3],
            ...     panel_labels=["Discovery", "Replication", "Meta-analysis"],
            ... )
        """
        n_gwas = len(gwas_dfs)
        if n_gwas == 0:
            raise ValueError("At least one GWAS DataFrame required")

        if panel_labels is not None and len(panel_labels) != n_gwas:
            raise ValueError(
                f"panel_labels length ({len(panel_labels)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )

        # Prepare all data first to get consistent x-axis
        prepared_dfs = []
        for df in gwas_dfs:
            prepared_df = prepare_manhattan_data(
                df=df,
                chrom_col=chrom_col,
                pos_col=pos_col,
                p_col=p_col,
                species=self.species,
                custom_order=custom_chrom_order,
            )
            prepared_dfs.append(prepared_df)

        # Use first df for chromosome order and centers
        chrom_order = prepared_dfs[0].attrs["chrom_order"]
        chrom_centers = prepared_dfs[0].attrs["chrom_centers"]

        # Calculate figure layout
        panel_height = figsize[1] / n_gwas
        height_ratios = [panel_height] * n_gwas

        # Create figure
        fig, axes = self._backend.create_figure(
            n_panels=n_gwas,
            height_ratios=height_ratios,
            figsize=figsize,
            sharex=True,
        )

        # Get consistent x limits across all panels
        x_min = min(df["_cumulative_pos"].min() for df in prepared_dfs)
        x_max = max(df["_cumulative_pos"].max() for df in prepared_dfs)
        x_padding = (x_max - x_min) * 0.01

        # Plot each panel
        for i, prepared_df in enumerate(prepared_dfs):
            ax = axes[i]

            # Plot points and significance line
            self._render_manhattan_points(ax, prepared_df, chrom_order)
            self._add_significance_line(ax, significance_threshold)

            # Set limits
            self._backend.set_xlim(ax, x_min - x_padding, x_max + x_padding)
            y_max = prepared_df["_neg_log_p"].max()
            self._backend.set_ylim(ax, 0, y_max * 1.1)

            # Labels
            self._backend.set_ylabel(ax, r"$-\log_{10}(p)$", fontsize=10)
            self._backend.hide_spines(ax, ["top", "right"])

            # Panel label
            if panel_labels and i < len(panel_labels):
                self._backend.add_panel_label(ax, panel_labels[i])

            # Set x-axis ticks for all panels (needed for interactive backends)
            positions = [
                chrom_centers[chrom] for chrom in chrom_order if chrom in chrom_centers
            ]
            labels = [chrom for chrom in chrom_order if chrom in chrom_centers]
            self._backend.set_xticks(ax, positions, labels, fontsize=8)

            # Only show x-axis label on bottom panel
            if i == n_gwas - 1:
                self._backend.set_xlabel(ax, "Chromosome", fontsize=12)

        # Overall title
        if title:
            self._backend.set_title(axes[0], title, fontsize=14)

        self._backend.finalize_layout(fig, hspace=0.1)

        return fig

    def plot_manhattan_qq(
        self,
        df: pd.DataFrame,
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        significance_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        show_confidence_band: bool = True,
        show_lambda: bool = True,
        figsize: Tuple[float, float] = (14, 5),
        title: Optional[str] = None,
    ) -> Any:
        """Create side-by-side Manhattan and QQ plots.

        Displays a Manhattan plot on the left and a QQ plot on the right,
        commonly used for GWAS publication figures.

        Args:
            df: GWAS results DataFrame.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            custom_chrom_order: Custom chromosome order (overrides species).
            significance_threshold: P-value threshold for genome-wide significance.
            show_confidence_band: If True, show 95% confidence band on QQ plot.
            show_lambda: If True, show genomic inflation factor on QQ plot.
            figsize: Figure size as (width, height).
            title: Overall plot title.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_manhattan_qq(gwas_df)
            >>> fig.savefig("gwas_summary.png", dpi=150)
        """
        # Prepare Manhattan data
        manhattan_df = prepare_manhattan_data(
            df=df,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )

        # Prepare QQ data
        qq_df = prepare_qq_data(df, p_col=p_col)

        # Create figure with side-by-side layout (Manhattan wider than QQ)
        # Use width_ratios to make Manhattan ~70% and QQ ~30%
        fig, axes = self._backend.create_figure_grid(
            n_rows=1,
            n_cols=2,
            width_ratios=[2.5, 1],
            figsize=figsize,
        )
        manhattan_ax = axes[0]
        qq_ax = axes[1]

        # --- Manhattan plot ---
        chrom_order = manhattan_df.attrs["chrom_order"]
        chrom_centers = manhattan_df.attrs["chrom_centers"]

        self._render_manhattan_points(manhattan_ax, manhattan_df, chrom_order)
        self._add_significance_line(manhattan_ax, significance_threshold)

        x_min = manhattan_df["_cumulative_pos"].min()
        x_max = manhattan_df["_cumulative_pos"].max()
        x_padding = (x_max - x_min) * 0.01
        self._backend.set_xlim(manhattan_ax, x_min - x_padding, x_max + x_padding)

        y_max = manhattan_df["_neg_log_p"].max()
        self._backend.set_ylim(manhattan_ax, 0, y_max * 1.1)

        positions = [
            chrom_centers[chrom] for chrom in chrom_order if chrom in chrom_centers
        ]
        labels = [chrom for chrom in chrom_order if chrom in chrom_centers]
        self._backend.set_xticks(manhattan_ax, positions, labels, fontsize=8)

        self._backend.set_xlabel(manhattan_ax, "Chromosome", fontsize=12)
        self._backend.set_ylabel(manhattan_ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_title(manhattan_ax, "Manhattan Plot", fontsize=12)
        self._backend.hide_spines(manhattan_ax, ["top", "right"])

        # --- QQ plot ---
        self._render_qq_plot(qq_ax, qq_df, show_confidence_band)

        self._backend.set_xlabel(qq_ax, r"Expected $-\log_{10}(p)$", fontsize=12)
        self._backend.set_ylabel(qq_ax, r"Observed $-\log_{10}(p)$", fontsize=12)

        if show_lambda:
            lambda_gc = qq_df.attrs["lambda_gc"]
            qq_title = f"QQ Plot (λ = {lambda_gc:.3f})"
        else:
            qq_title = "QQ Plot"
        self._backend.set_title(qq_ax, qq_title, fontsize=12)
        self._backend.hide_spines(qq_ax, ["top", "right"])

        # Overall title - adjust top margin to leave room for suptitle
        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90)
        else:
            self._backend.finalize_layout(fig)

        return fig

    def plot_manhattan_qq_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        significance_threshold: Optional[float] = DEFAULT_GENOMEWIDE_THRESHOLD,
        show_confidence_band: bool = True,
        show_lambda: bool = True,
        panel_labels: Optional[List[str]] = None,
        figsize: Tuple[float, float] = (14, 8),
        title: Optional[str] = None,
    ) -> Any:
        """Create stacked side-by-side Manhattan and QQ plots for multiple GWAS.

        Displays Manhattan+QQ pairs for each GWAS dataset, stacked vertically
        for easy comparison across studies.

        Args:
            gwas_dfs: List of GWAS results DataFrames.
            chrom_col: Column name for chromosome.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            custom_chrom_order: Custom chromosome order (overrides species).
            significance_threshold: P-value threshold for genome-wide significance.
            show_confidence_band: If True, show 95% confidence band on QQ plots.
            show_lambda: If True, show genomic inflation factor on QQ plots.
            panel_labels: List of labels for each GWAS (one per dataset).
            figsize: Figure size as (width, height).
            title: Overall plot title.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_manhattan_qq_stacked(
            ...     [discovery_df, replication_df],
            ...     panel_labels=["Discovery", "Replication"],
            ... )
        """
        n_gwas = len(gwas_dfs)
        if n_gwas == 0:
            raise ValueError("At least one GWAS DataFrame required")

        # Prepare all data
        manhattan_dfs = []
        qq_dfs = []
        for df in gwas_dfs:
            manhattan_dfs.append(
                prepare_manhattan_data(
                    df=df,
                    chrom_col=chrom_col,
                    pos_col=pos_col,
                    p_col=p_col,
                    species=self.species,
                    custom_order=custom_chrom_order,
                )
            )
            qq_dfs.append(prepare_qq_data(df, p_col=p_col))

        # Use chromosome order from first dataset
        chrom_order = manhattan_dfs[0].attrs["chrom_order"]
        chrom_centers = manhattan_dfs[0].attrs["chrom_centers"]

        # Create grid: n_gwas rows, 2 columns (Manhattan | QQ)
        # Manhattan gets ~70% width, QQ gets ~30%
        fig, axes = self._backend.create_figure_grid(
            n_rows=n_gwas,
            n_cols=2,
            width_ratios=[2.5, 1],
            figsize=figsize,
        )

        # Get consistent x limits for Manhattan plots
        x_min = min(df["_cumulative_pos"].min() for df in manhattan_dfs)
        x_max = max(df["_cumulative_pos"].max() for df in manhattan_dfs)
        x_padding = (x_max - x_min) * 0.01

        # Plot each row
        for i in range(n_gwas):
            manhattan_ax = axes[i * 2]  # Even indices: Manhattan
            qq_ax = axes[i * 2 + 1]  # Odd indices: QQ
            manhattan_df = manhattan_dfs[i]
            qq_df = qq_dfs[i]

            # --- Manhattan plot ---
            self._render_manhattan_points(manhattan_ax, manhattan_df, chrom_order)
            self._add_significance_line(manhattan_ax, significance_threshold)

            self._backend.set_xlim(manhattan_ax, x_min - x_padding, x_max + x_padding)
            y_max = manhattan_df["_neg_log_p"].max()
            self._backend.set_ylim(manhattan_ax, 0, y_max * 1.1)

            # Panel label
            if panel_labels and i < len(panel_labels):
                self._backend.add_panel_label(manhattan_ax, panel_labels[i])

            # Y-axis label
            self._backend.set_ylabel(manhattan_ax, r"$-\log_{10}(p)$", fontsize=10)
            self._backend.hide_spines(manhattan_ax, ["top", "right"])

            # X-axis: set chromosome ticks for all panels (consistent with plot_manhattan_stacked)
            positions = [
                chrom_centers[chrom] for chrom in chrom_order if chrom in chrom_centers
            ]
            chrom_labels = [chrom for chrom in chrom_order if chrom in chrom_centers]
            self._backend.set_xticks(manhattan_ax, positions, chrom_labels, fontsize=8)

            # Only show "Chromosome" label on bottom row
            if i == n_gwas - 1:
                self._backend.set_xlabel(manhattan_ax, "Chromosome", fontsize=10)

            # --- QQ plot ---
            self._render_qq_plot(qq_ax, qq_df, show_confidence_band)

            # Labels for QQ
            if i == n_gwas - 1:
                self._backend.set_xlabel(
                    qq_ax, r"Expected $-\log_{10}(p)$", fontsize=10
                )
            self._backend.set_ylabel(qq_ax, r"Observed $-\log_{10}(p)$", fontsize=10)

            # QQ title with lambda
            if show_lambda:
                lambda_gc = qq_df.attrs["lambda_gc"]
                qq_title = f"λ = {lambda_gc:.3f}"
            else:
                qq_title = "QQ"
            self._backend.set_title(qq_ax, qq_title, fontsize=10)
            self._backend.hide_spines(qq_ax, ["top", "right"])

        # Overall title - adjust top margin to leave room for suptitle
        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90, hspace=0.15)
        else:
            self._backend.finalize_layout(fig, hspace=0.15)

        return fig

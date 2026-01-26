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
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator

from .backends import BackendType, get_backend
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
)
from .eqtl import validate_eqtl_df
from .finemapping import (
    get_credible_sets,
    prepare_finemapping_for_plotting,
)
from .gene_track import assign_gene_positions, plot_gene_track
from .labels import add_snp_labels
from .ld import calculate_ld, find_plink
from .logging import enable_logging, logger
from .recombination import (
    add_recombination_overlay,
    download_canine_recombination_maps,
    get_default_data_dir,
    get_recombination_rate_for_region,
)
from .utils import normalize_chrom, validate_genes_df, validate_gwas_df

# Default significance threshold: 5e-8 for human, 5e-7 for canine
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-7
DEFAULT_GENOMEWIDE_LINE = -np.log10(DEFAULT_GENOMEWIDE_THRESHOLD)


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
    ):
        """Initialize the plotter."""
        # Configure logging
        if log_level is not None:
            enable_logging(log_level)

        self.species = species
        self.genome_build = (
            genome_build if genome_build else self._default_build(species)
        )
        self.backend_name = backend
        self._backend = get_backend(backend)
        self.plink_path = plink_path or find_plink()
        self.recomb_data_dir = recomb_data_dir
        self.genomewide_threshold = genomewide_threshold
        self._genomewide_line = -np.log10(genomewide_threshold)

        # Cache for loaded data
        self._recomb_cache = {}

    @staticmethod
    def _default_build(species: str) -> Optional[str]:
        """Get default genome build for species."""
        if species == "canine":
            return "canfam3.1"
        if species == "feline":
            return "felCat9"
        return None

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
                and len(list(default_dir.glob("chr*_recomb.tsv"))) >= 38
            ):
                return default_dir
            # Download
            try:
                return download_canine_recombination_maps()
            except Exception as e:
                logger.warning(f"Could not download recombination maps: {e}")
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

    def plot(
        self,
        gwas_df: pd.DataFrame,
        chrom: int,
        start: int,
        end: int,
        lead_pos: Optional[int] = None,
        ld_reference_file: Optional[str] = None,
        ld_col: Optional[str] = None,
        genes_df: Optional[pd.DataFrame] = None,
        exons_df: Optional[pd.DataFrame] = None,
        recomb_df: Optional[pd.DataFrame] = None,
        show_recombination: bool = True,
        snp_labels: bool = True,
        label_top_n: int = 5,
        pos_col: str = "ps",
        p_col: str = "p_wald",
        rs_col: str = "rs",
        figsize: Tuple[int, int] = (12, 8),
    ) -> Figure:
        """Create a regional association plot.

        Args:
            gwas_df: GWAS results DataFrame.
            chrom: Chromosome number.
            start: Start position of the region.
            end: End position of the region.
            lead_pos: Position of the lead/index SNP to highlight.
            ld_reference_file: PLINK binary fileset for LD calculation.
                If provided with lead_pos, calculates LD on the fly.
            ld_col: Column name for pre-computed LD (R²) values.
                Use this if LD was calculated externally.
            genes_df: Gene annotations with chr, start, end, gene_name.
            exons_df: Exon annotations with chr, start, end, gene_name.
            recomb_df: Pre-loaded recombination rate data.
                If None and show_recombination=True, loads from species default.
            show_recombination: Whether to show recombination rate overlay.
            snp_labels: Whether to label top SNPs.
            label_top_n: Number of top SNPs to label.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP ID.
            figsize: Figure size.

        Returns:
            Matplotlib Figure object.

        Raises:
            ValidationError: If required DataFrame columns are missing.
        """
        # Validate inputs
        validate_gwas_df(gwas_df, pos_col=pos_col, p_col=p_col)
        if genes_df is not None:
            validate_genes_df(genes_df)

        logger.debug(f"Creating plot for chr{chrom}:{start}-{end}")

        # Prevent auto-display in interactive environments
        plt.ioff()

        # Prepare data
        df = gwas_df.copy()
        df["neglog10p"] = -np.log10(df[p_col].clip(lower=1e-300))

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
        self._plot_association(ax, df, pos_col, ld_col, lead_pos)

        # Add significance line
        ax.axhline(
            y=self._genomewide_line,
            color="red",
            linestyle=(0, (5, 10)),
            linewidth=1,
            alpha=0.8,
            zorder=1,
        )

        # Add SNP labels
        if snp_labels and rs_col in df.columns and label_top_n > 0 and not df.empty:
            add_snp_labels(
                ax,
                df,
                pos_col=pos_col,
                neglog10p_col="neglog10p",
                rs_col=rs_col,
                label_top_n=label_top_n,
                genes_df=genes_df,
                chrom=chrom,
            )

        # Add recombination overlay
        if recomb_df is not None and not recomb_df.empty:
            add_recombination_overlay(ax, recomb_df, start, end)

        # Format axes
        ax.set_ylabel(r"$-\log_{10}$ P")
        ax.set_xlim(start, end)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Add LD legend
        if ld_col is not None and ld_col in df.columns:
            self._add_ld_legend(ax)

        # Plot gene track
        if genes_df is not None and gene_ax is not None:
            plot_gene_track(gene_ax, genes_df, chrom, start, end, exons_df)
            gene_ax.set_xlabel(f"Chromosome {chrom} (Mb)")
            gene_ax.spines["top"].set_visible(False)
            gene_ax.spines["right"].set_visible(False)
            gene_ax.spines["left"].set_visible(False)
        else:
            ax.set_xlabel(f"Chromosome {chrom} (Mb)")

        # Format x-axis with Mb labels
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:.2f}"))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))

        # Adjust layout
        fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.1, hspace=0.08)
        plt.ion()

        return fig

    def _create_figure(
        self,
        genes_df: Optional[pd.DataFrame],
        chrom: int,
        start: int,
        end: int,
        figsize: Tuple[int, int],
    ) -> Tuple[Figure, Axes, Optional[Axes]]:
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

            fig, axes = plt.subplots(
                2,
                1,
                figsize=(figsize[0], total_height),
                height_ratios=[assoc_height, gene_track_height],
                sharex=True,
                gridspec_kw={"hspace": 0},
            )
            return fig, axes[0], axes[1]
        else:
            fig, ax = plt.subplots(figsize=(figsize[0], figsize[1] * 0.75))
            return fig, ax, None

    def _plot_association(
        self,
        ax: Axes,
        df: pd.DataFrame,
        pos_col: str,
        ld_col: Optional[str],
        lead_pos: Optional[int],
    ) -> None:
        """Plot association scatter with LD coloring."""
        # LD-based coloring
        if ld_col is not None and ld_col in df.columns:
            df["ld_bin"] = df[ld_col].apply(get_ld_bin)
            df = df.sort_values(ld_col, ascending=True, na_position="first")

            palette = get_ld_color_palette()
            for bin_label in df["ld_bin"].unique():
                bin_data = df[df["ld_bin"] == bin_label]
                ax.scatter(
                    bin_data[pos_col],
                    bin_data["neglog10p"],
                    c=palette.get(bin_label, "#BEBEBE"),
                    s=60,
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                )
        else:
            # Default: grey points
            ax.scatter(
                df[pos_col],
                df["neglog10p"],
                c="#BEBEBE",
                s=60,
                edgecolor="black",
                linewidth=0.5,
                zorder=2,
            )

        # Highlight lead SNP
        if lead_pos is not None:
            lead_snp = df[df[pos_col] == lead_pos]
            if not lead_snp.empty:
                ax.scatter(
                    lead_snp[pos_col],
                    lead_snp["neglog10p"],
                    c=LEAD_SNP_COLOR,
                    s=60,
                    marker="D",
                    edgecolors="black",
                    linewidths=1.5,
                    zorder=10,
                )

    def _add_ld_legend(self, ax: Axes) -> None:
        """Add LD color legend to plot."""
        palette = get_ld_color_palette()
        legend_elements = [
            Line2D(
                [0],
                [0],
                marker="D",
                color="w",
                markerfacecolor=LEAD_SNP_COLOR,
                markeredgecolor="black",
                markersize=6,
                label="Lead SNP",
            ),
        ]

        for threshold, label, _ in LD_BINS:
            legend_elements.append(
                Patch(
                    facecolor=palette[label],
                    edgecolor="black",
                    label=label,
                )
            )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=9,
            frameon=True,
            framealpha=0.9,
            title=r"$r^2$",
            title_fontsize=10,
            handlelength=1.5,
            handleheight=1.0,
            labelspacing=0.4,
        )

    def _add_eqtl_legend(self, ax: Axes) -> None:
        """Add eQTL effect size legend to plot."""
        legend_elements = []

        # Positive effects (upward triangles)
        for _, _, label, color in EQTL_POSITIVE_BINS:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="^",
                    color="w",
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markersize=7,
                    label=label,
                )
            )

        # Negative effects (downward triangles)
        for _, _, label, color in EQTL_NEGATIVE_BINS:
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="v",
                    color="w",
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markersize=7,
                    label=label,
                )
            )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=8,
            frameon=True,
            framealpha=0.9,
            title="eQTL effect",
            title_fontsize=9,
            handlelength=1.2,
            handleheight=1.0,
            labelspacing=0.3,
        )

    def _plot_finemapping(
        self,
        ax: Axes,
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
        # Sort by position for line plotting
        df = df.sort_values(pos_col)

        # Plot PIP as line
        ax.plot(
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
                ax.scatter(
                    cs_data[pos_col],
                    cs_data[pip_col],
                    c=color,
                    s=50,
                    marker="o",
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=3,
                    label=f"CS{cs_id}",
                )
            # Plot variants not in any credible set
            non_cs_data = df[(df[cs_col].isna()) | (df[cs_col] == 0)]
            if not non_cs_data.empty and pip_threshold > 0:
                non_cs_data = non_cs_data[non_cs_data[pip_col] >= pip_threshold]
                if not non_cs_data.empty:
                    ax.scatter(
                        non_cs_data[pos_col],
                        non_cs_data[pip_col],
                        c="#BEBEBE",
                        s=30,
                        marker="o",
                        edgecolor="black",
                        linewidth=0.3,
                        zorder=2,
                        alpha=0.6,
                    )
        else:
            # No credible sets - show all points above threshold
            if pip_threshold > 0:
                high_pip = df[df[pip_col] >= pip_threshold]
                if not high_pip.empty:
                    ax.scatter(
                        high_pip[pos_col],
                        high_pip[pip_col],
                        c=PIP_LINE_COLOR,
                        s=50,
                        marker="o",
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=3,
                    )

    def _add_finemapping_legend(
        self,
        ax: Axes,
        credible_sets: List[int],
    ) -> None:
        """Add fine-mapping legend showing credible sets.

        Args:
            ax: Matplotlib axes object.
            credible_sets: List of credible set IDs to include.
        """
        if not credible_sets:
            return

        legend_elements = []
        for cs_id in credible_sets:
            color = get_credible_set_color(cs_id)
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=color,
                    markeredgecolor="black",
                    markersize=7,
                    label=f"CS{cs_id}",
                )
            )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            fontsize=8,
            frameon=True,
            framealpha=0.9,
            title="Credible sets",
            title_fontsize=9,
            handlelength=1.2,
            handleheight=1.0,
            labelspacing=0.3,
        )

    def plot_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        chrom: int,
        start: int,
        end: int,
        lead_positions: Optional[List[int]] = None,
        panel_labels: Optional[List[str]] = None,
        ld_reference_file: Optional[str] = None,
        ld_reference_files: Optional[List[str]] = None,
        ld_col: Optional[str] = None,
        genes_df: Optional[pd.DataFrame] = None,
        exons_df: Optional[pd.DataFrame] = None,
        eqtl_df: Optional[pd.DataFrame] = None,
        eqtl_gene: Optional[str] = None,
        finemapping_df: Optional[pd.DataFrame] = None,
        finemapping_cs_col: Optional[str] = "cs",
        recomb_df: Optional[pd.DataFrame] = None,
        show_recombination: bool = True,
        snp_labels: bool = True,
        label_top_n: int = 3,
        pos_col: str = "ps",
        p_col: str = "p_wald",
        rs_col: str = "rs",
        figsize: Tuple[float, Optional[float]] = (12, None),
    ) -> Any:
        """Create stacked regional association plots for multiple GWAS.

        Vertically stacks multiple GWAS results for comparison, with shared
        x-axis and optional gene track at the bottom.

        Args:
            gwas_dfs: List of GWAS results DataFrames to stack.
            chrom: Chromosome number.
            start: Start position of the region.
            end: End position of the region.
            lead_positions: List of lead SNP positions (one per GWAS).
                If None, auto-detects from lowest p-value.
            panel_labels: Labels for each panel (e.g., phenotype names).
            ld_reference_file: Single PLINK fileset for all panels.
            ld_reference_files: List of PLINK filesets (one per panel).
            ld_col: Column name for pre-computed LD (R²) values in each DataFrame.
                Use this if LD was calculated externally.
            genes_df: Gene annotations for bottom track.
            exons_df: Exon annotations for gene track.
            eqtl_df: eQTL data to display as additional panel.
            eqtl_gene: Filter eQTL data to this target gene.
            finemapping_df: Fine-mapping/SuSiE results with pos and pip columns.
                Displayed as PIP line with optional credible set coloring.
            finemapping_cs_col: Column name for credible set assignment in finemapping_df.
            recomb_df: Pre-loaded recombination rate data.
            show_recombination: Whether to show recombination overlay.
            snp_labels: Whether to label top SNPs.
            label_top_n: Number of top SNPs to label per panel.
            pos_col: Column name for position.
            p_col: Column name for p-value.
            rs_col: Column name for SNP ID.
            figsize: Figure size (width, height). If height is None, auto-calculates.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_stacked(
            ...     [gwas_height, gwas_bmi, gwas_whr],
            ...     chrom=1, start=1000000, end=2000000,
            ...     panel_labels=["Height", "BMI", "WHR"],
            ...     genes_df=genes_df,
            ... )
        """
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
                    lead_idx = region_df[p_col].idxmin()
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

        # Prevent auto-display in interactive environments
        plt.ioff()

        # Load recombination data if needed
        if show_recombination and recomb_df is None:
            recomb_df = self._get_recomb_for_region(chrom, start, end)

        # Create figure
        fig, axes = plt.subplots(
            n_panels,
            1,
            figsize=actual_figsize,
            height_ratios=height_ratios,
            sharex=True,
            gridspec_kw={"hspace": 0.05},
        )
        if n_panels == 1:
            axes = [axes]

        # Plot each GWAS panel
        for i, (gwas_df, lead_pos) in enumerate(zip(gwas_dfs, lead_positions)):
            ax = axes[i]
            df = gwas_df.copy()
            df["neglog10p"] = -np.log10(df[p_col].clip(lower=1e-300))

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
            self._plot_association(ax, df, pos_col, panel_ld_col, lead_pos)

            # Add significance line
            ax.axhline(
                y=self._genomewide_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.8,
                zorder=1,
            )

            # Add SNP labels
            if snp_labels and rs_col in df.columns and label_top_n > 0 and not df.empty:
                add_snp_labels(
                    ax,
                    df,
                    pos_col=pos_col,
                    neglog10p_col="neglog10p",
                    rs_col=rs_col,
                    label_top_n=label_top_n,
                    genes_df=genes_df,
                    chrom=chrom,
                )

            # Add recombination overlay (only on first panel)
            if i == 0 and recomb_df is not None and not recomb_df.empty:
                add_recombination_overlay(ax, recomb_df, start, end)

            # Format axes
            ax.set_ylabel(r"$-\log_{10}$ P")
            ax.set_xlim(start, end)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # Add panel label
            if panel_labels and i < len(panel_labels):
                ax.annotate(
                    panel_labels[i],
                    xy=(0.02, 0.95),
                    xycoords="axes fraction",
                    fontsize=11,
                    fontweight="bold",
                    va="top",
                    ha="left",
                )

            # Add LD legend (only on first panel)
            if i == 0 and panel_ld_col is not None and panel_ld_col in df.columns:
                self._add_ld_legend(ax)

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

                # Add legend for credible sets
                credible_sets = get_credible_sets(fm_data, finemapping_cs_col)
                if credible_sets:
                    self._add_finemapping_legend(ax, credible_sets)

            ax.set_ylabel("PIP")
            ax.set_ylim(-0.05, 1.05)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            panel_idx += 1

        # Plot eQTL panel if provided
        eqtl_panel_idx = panel_idx
        if eqtl_df is not None:
            ax = axes[eqtl_panel_idx]
            eqtl_data = eqtl_df.copy()

            # Filter by gene if specified
            if eqtl_gene and "gene" in eqtl_data.columns:
                eqtl_data = eqtl_data[eqtl_data["gene"] == eqtl_gene]

            # Filter by region
            if "pos" in eqtl_data.columns:
                eqtl_data = eqtl_data[
                    (eqtl_data["pos"] >= start) & (eqtl_data["pos"] <= end)
                ]

            if not eqtl_data.empty:
                eqtl_data["neglog10p"] = -np.log10(
                    eqtl_data["p_value"].clip(lower=1e-300)
                )

                # Check if effect_size column exists for directional coloring
                has_effect = "effect_size" in eqtl_data.columns

                if has_effect:
                    # Plot triangles by effect direction with color by magnitude
                    for _, row in eqtl_data.iterrows():
                        effect = row["effect_size"]
                        color = get_eqtl_color(effect)
                        marker = "^" if effect >= 0 else "v"
                        ax.scatter(
                            row["pos"],
                            row["neglog10p"],
                            c=color,
                            s=50,
                            marker=marker,
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
                    # Add eQTL effect legend
                    self._add_eqtl_legend(ax)
                else:
                    # No effect sizes - plot as diamonds
                    ax.scatter(
                        eqtl_data["pos"],
                        eqtl_data["neglog10p"],
                        c="#FF6B6B",
                        s=60,
                        marker="D",
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=2,
                        label=f"eQTL ({eqtl_gene})" if eqtl_gene else "eQTL",
                    )
                    ax.legend(loc="upper right", fontsize=9)

            ax.set_ylabel(r"$-\log_{10}$ P (eQTL)")
            ax.axhline(
                y=self._genomewide_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.8,
            )
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            panel_idx += 1

        # Plot gene track
        if genes_df is not None:
            gene_ax = axes[panel_idx]
            plot_gene_track(gene_ax, genes_df, chrom, start, end, exons_df)
            gene_ax.set_xlabel(f"Chromosome {chrom} (Mb)")
            gene_ax.spines["top"].set_visible(False)
            gene_ax.spines["right"].set_visible(False)
            gene_ax.spines["left"].set_visible(False)
        else:
            # Set x-label on bottom panel
            axes[-1].set_xlabel(f"Chromosome {chrom} (Mb)")

        # Format x-axis
        axes[0].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:.2f}"))
        axes[0].xaxis.set_major_locator(MaxNLocator(nbins=6))

        # Adjust layout
        fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08, hspace=0.05)
        plt.ion()

        return fig

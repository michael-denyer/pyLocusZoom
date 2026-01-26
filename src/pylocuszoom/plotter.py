"""Main LocusZoomPlotter class for regional association plots.

Orchestrates all components (LD coloring, gene track, recombination overlay,
SNP labels) into a unified plotting interface.
"""

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator

from .colors import (
    LD_BINS,
    LEAD_SNP_COLOR,
    get_ld_bin,
    get_ld_color_palette,
)
from .gene_track import assign_gene_positions, plot_gene_track
from .labels import add_snp_labels
from .ld import calculate_ld, find_plink
from .logging import enable_logging, logger
from .recombination import (
    add_recombination_overlay,
    download_dog_recombination_maps,
    get_default_data_dir,
    get_recombination_rate_for_region,
)
from .utils import normalize_chrom, validate_genes_df, validate_gwas_df

# Default significance threshold: 5e-8 for human, 5e-7 for dog
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-7
DEFAULT_GENOMEWIDE_LINE = -np.log10(DEFAULT_GENOMEWIDE_THRESHOLD)


class LocusZoomPlotter:
    """Regional association plot generator with LD coloring and annotations.

    Creates LocusZoom-style regional plots with:
    - LD coloring based on R² with lead variant
    - Gene and exon tracks
    - Recombination rate overlays (dog built-in, or user-provided)
    - Automatic SNP labeling

    Args:
        species: Species name ('dog', 'cat', or None for custom).
            Dog has built-in recombination maps.
        genome_build: Genome build for coordinate system. For dog:
            "canfam3.1" (default) or "canfam4". If "canfam4", recombination
            maps are automatically lifted over from CanFam3.1.
        plink_path: Path to PLINK executable for LD calculation.
            Auto-detects if None.
        recomb_data_dir: Directory containing recombination maps.
            Uses platform cache if None.
        genomewide_threshold: P-value threshold for significance line.
        log_level: Logging level ("DEBUG", "INFO", "WARNING", "ERROR", or None
            to disable). Defaults to "INFO".

    Example:
        >>> # CanFam3.1 (default)
        >>> plotter = LocusZoomPlotter(species="dog")
        >>>
        >>> # CanFam4
        >>> plotter = LocusZoomPlotter(species="dog", genome_build="canfam4")
        >>>
        >>> fig = plotter.plot(
        ...     gwas_df,
        ...     chrom=1,
        ...     start=1000000,
        ...     end=2000000,
        ...     lead_pos=1500000,
        ... )
        >>> fig.savefig("regional_plot.png", dpi=150)
    """

    def __init__(
        self,
        species: str = "dog",
        genome_build: Optional[str] = None,
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
        self.plink_path = plink_path or find_plink()
        self.recomb_data_dir = recomb_data_dir
        self.genomewide_threshold = genomewide_threshold
        self._genomewide_line = -np.log10(genomewide_threshold)

        # Cache for loaded data
        self._recomb_cache = {}

    @staticmethod
    def _default_build(species: str) -> Optional[str]:
        """Get default genome build for species."""
        if species == "dog":
            return "canfam3.1"
        if species == "cat":
            return "felCat9"
        return None

    def _ensure_recomb_maps(self) -> Optional[Path]:
        """Ensure recombination maps are downloaded.

        Returns path to recombination map directory, or None if not available.
        """
        if self.species == "dog":
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
                return download_dog_recombination_maps()
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
            color="grey",
            linestyle="--",
            linewidth=1,
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
                    s=120,
                    marker="D",
                    edgecolors="black",
                    linewidths=1,
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
                markersize=8,
                label="Index SNP",
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
            loc="upper left",
            fontsize=9,
            frameon=True,
            framealpha=0.9,
            title=r"$r^2$",
            title_fontsize=10,
            handlelength=1.5,
            handleheight=1.0,
            labelspacing=0.4,
        )

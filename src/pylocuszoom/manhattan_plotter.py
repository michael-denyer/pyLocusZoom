"""Manhattan and QQ plot generator.

Provides genome-wide visualization of GWAS results including:
- Manhattan plots (standard and categorical)
- QQ plots with confidence bands
- Combined Manhattan+QQ layouts
- Stacked multi-GWAS comparisons
"""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._plotter_utils import add_significance_line
from .backends import BackendType, get_backend
from .manhattan import prepare_categorical_data, prepare_manhattan_data
from .qq import prepare_qq_data

# Manhattan/QQ plot styling constants
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8
MANHATTAN_POINT_SIZE = 10
MANHATTAN_CATEGORICAL_POINT_SIZE = 30
QQ_POINT_SIZE = 10
POINT_EDGE_COLOR = "black"
MANHATTAN_EDGE_WIDTH = 0.1
QQ_EDGE_WIDTH = 0.02
QQ_POINT_COLOR = "#1f77b4"
QQ_CI_COLOR = "#CCCCCC"
QQ_CI_ALPHA = 0.5
SIGNIFICANCE_LINE_COLOR = "red"


class ManhattanPlotter:
    """Manhattan and QQ plot generator for genome-wide visualizations.

    Creates publication-quality Manhattan plots, QQ plots, and combined
    layouts for GWAS summary statistics.

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        species: Species name ('canine', 'feline', 'human', or None).
            Used to determine chromosome order.
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').
        genomewide_threshold: P-value threshold for significance line.

    Example:
        >>> plotter = ManhattanPlotter(species="human")
        >>> fig = plotter.plot_manhattan(gwas_df)
        >>> fig.savefig("manhattan.png", dpi=150)
    """

    def __init__(
        self,
        species: str = "canine",
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
    ):
        """Initialize the Manhattan plotter."""
        self.species = species
        self._backend = get_backend(backend)
        self.genomewide_threshold = genomewide_threshold

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
        add_significance_line(self._backend, ax, significance_threshold)

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

        add_significance_line(self._backend, ax, significance_threshold)

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
            add_significance_line(self._backend, ax, significance_threshold)

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
        add_significance_line(self._backend, manhattan_ax, significance_threshold)

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

        # Overall title
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
            add_significance_line(self._backend, manhattan_ax, significance_threshold)

            self._backend.set_xlim(manhattan_ax, x_min - x_padding, x_max + x_padding)
            y_max = manhattan_df["_neg_log_p"].max()
            self._backend.set_ylim(manhattan_ax, 0, y_max * 1.1)

            # Panel label
            if panel_labels and i < len(panel_labels):
                self._backend.add_panel_label(manhattan_ax, panel_labels[i])

            # Y-axis label
            self._backend.set_ylabel(manhattan_ax, r"$-\log_{10}(p)$", fontsize=10)
            self._backend.hide_spines(manhattan_ax, ["top", "right"])

            # X-axis: set chromosome ticks for all panels
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

        # Overall title
        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90, hspace=0.15)
        else:
            self._backend.finalize_layout(fig, hspace=0.15)

        return fig

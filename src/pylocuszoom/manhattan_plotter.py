"""Manhattan and QQ plot generator.

Provides genome-wide visualization of GWAS results including:
- Manhattan plots (standard and categorical)
- QQ plots with confidence bands
- Combined Manhattan+QQ layouts
- Stacked multi-GWAS comparisons
"""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from ._rendering import ManhattanQQRenderer
from .backends import BackendType, get_backend
from .manhattan import (
    prepare_categorical_data,
    prepare_manhattan_data,
    prepare_manhattan_frames,
)
from .qq import prepare_qq_data
from .species import Species, resolve_species


class ManhattanPlotter:
    """Manhattan and QQ plot generator for genome-wide visualizations.

    Creates publication-quality Manhattan plots, QQ plots, and combined
    layouts for GWAS summary statistics.

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        species: Species name, alias or record ('canine', 'dog', 'feline',
            'human', or None). An unknown name raises ValidationError.
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
        species: str | Species | None = "canine",
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
    ):
        """Initialize the Manhattan plotter."""
        self.species = resolve_species(species)
        self._backend = get_backend(backend)
        self._renderer = ManhattanQQRenderer(self._backend)
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
        significance_threshold: ThresholdArg = UNSET,
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
            significance_threshold: P-value threshold for the genome-wide
                significance line. Defaults to the plotter's
                ``genomewide_threshold``; pass None to draw no line.
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
        significance_threshold = resolve_threshold(
            significance_threshold, self.genomewide_threshold
        )

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

        return self._renderer.render_manhattan(
            prepared_df,
            figsize=figsize,
            significance_threshold=significance_threshold,
            title=title,
        )

    def _plot_manhattan_categorical(
        self,
        df: pd.DataFrame,
        category_col: str,
        p_col: str = "p",
        category_order: Optional[List[str]] = None,
        significance_threshold: Optional[float] = None,
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
        return self._renderer.render_categorical(
            prepared_df,
            figsize=figsize,
            significance_threshold=significance_threshold,
            title=title,
        )

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
        return self._renderer.render_qq(
            prepared_df,
            figsize=figsize,
            show_confidence_band=show_confidence_band,
            show_lambda=show_lambda,
            title=title,
        )

    def plot_manhattan_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        significance_threshold: ThresholdArg = UNSET,
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
            significance_threshold: P-value threshold for the genome-wide
                significance line. Defaults to the plotter's
                ``genomewide_threshold``; pass None to draw no line.
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
        significance_threshold = resolve_threshold(
            significance_threshold, self.genomewide_threshold
        )
        n_gwas = len(gwas_dfs)
        if n_gwas == 0:
            raise ValueError("At least one GWAS DataFrame required")

        if panel_labels is not None and len(panel_labels) != n_gwas:
            raise ValueError(
                f"panel_labels length ({len(panel_labels)}) must match "
                f"number of GWAS DataFrames ({n_gwas})"
            )

        prepared_dfs = prepare_manhattan_frames(
            gwas_dfs,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )
        return self._renderer.render_manhattan_stacked(
            prepared_dfs,
            figsize=figsize,
            significance_threshold=significance_threshold,
            panel_labels=panel_labels,
            title=title,
        )

    def plot_manhattan_qq(
        self,
        df: pd.DataFrame,
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        significance_threshold: ThresholdArg = UNSET,
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
            significance_threshold: P-value threshold for the genome-wide
                significance line. Defaults to the plotter's
                ``genomewide_threshold``; pass None to draw no line.
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
        significance_threshold = resolve_threshold(
            significance_threshold, self.genomewide_threshold
        )

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
        return self._renderer.render_manhattan_qq(
            manhattan_df,
            qq_df,
            figsize=figsize,
            significance_threshold=significance_threshold,
            show_confidence_band=show_confidence_band,
            show_lambda=show_lambda,
            title=title,
        )

    def plot_manhattan_qq_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        chrom_col: str = "chrom",
        pos_col: str = "pos",
        p_col: str = "p",
        custom_chrom_order: Optional[List[str]] = None,
        significance_threshold: ThresholdArg = UNSET,
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
            significance_threshold: P-value threshold for the genome-wide
                significance line. Defaults to the plotter's
                ``genomewide_threshold``; pass None to draw no line.
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
        significance_threshold = resolve_threshold(
            significance_threshold, self.genomewide_threshold
        )
        n_gwas = len(gwas_dfs)
        if n_gwas == 0:
            raise ValueError("At least one GWAS DataFrame required")

        manhattan_dfs = prepare_manhattan_frames(
            gwas_dfs,
            chrom_col=chrom_col,
            pos_col=pos_col,
            p_col=p_col,
            species=self.species,
            custom_order=custom_chrom_order,
        )
        qq_dfs = [prepare_qq_data(df, p_col=p_col) for df in gwas_dfs]
        return self._renderer.render_manhattan_qq_stacked(
            manhattan_dfs,
            qq_dfs,
            figsize=figsize,
            significance_threshold=significance_threshold,
            show_confidence_band=show_confidence_band,
            show_lambda=show_lambda,
            panel_labels=panel_labels,
            title=title,
        )

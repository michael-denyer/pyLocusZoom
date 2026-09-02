"""Manhattan and QQ plot generator.

Provides genome-wide visualization of GWAS results including:
- Manhattan plots (standard and categorical)
- QQ plots with confidence bands
- Combined Manhattan+QQ layouts
- Stacked multi-GWAS comparisons
"""

from typing import Any, List, Optional, Tuple

import pandas as pd

from ._figure import FigurePlan, render_figure
from ._manhattan_panel import (
    categorical_spec,
    manhattan_spec,
    stacked_manhattan_specs,
)
from ._plotter_utils import (
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from ._qq_panel import QQPanelSpec, qq_title
from .backends import BackendType, get_backend
from .config import GenomeWideConfig
from .manhattan import prepare_categorical_data, prepare_genomewide_frames
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

    Every method takes a :class:`~pylocuszoom.GenomeWideConfig` naming the
    chromosome, position and p-value columns and an optional chromosome
    order, and validates each frame against it before laying it out.

    Example:
        >>> plotter = ManhattanPlotter(species="human")
        >>> fig = plotter.plot_manhattan(gwas_df)
        >>> fig.savefig("manhattan.png", dpi=150)
        >>>
        >>> from pylocuszoom import GenomeWideConfig
        >>> columns = GenomeWideConfig(chrom_col="CHR", pos_col="BP", p_col="P")
        >>> fig = plotter.plot_manhattan(plink_df, config=columns)
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
        self.genomewide_threshold = genomewide_threshold

    def plot_manhattan(
        self,
        df: pd.DataFrame,
        *,
        config: GenomeWideConfig = GenomeWideConfig(),
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
            config: Column names and chromosome order. A categorical plot
                reads only ``config.p_col``.
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

        Raises:
            ValidationError: If ``df`` is empty or lacks a configured column.

        Example:
            >>> # Standard Manhattan plot
            >>> fig = plotter.plot_manhattan(gwas_df)
            >>>
            >>> # Categorical Manhattan (PheWAS-style)
            >>> fig = plotter.plot_manhattan(
            ...     phewas_df,
            ...     category_col="phenotype_category",
            ...     config=GenomeWideConfig(p_col="pvalue"),
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
                p_col=config.p_col,
                category_order=category_order,
                significance_threshold=significance_threshold,
                figsize=figsize,
                title=title,
            )

        # Standard Manhattan plot
        prepared = prepare_genomewide_frames([df], config, species=self.species)[0]

        panel = manhattan_spec(
            prepared,
            significance_threshold=significance_threshold,
            x_label="Chromosome",
            title=title or "Manhattan Plot",
        )
        return render_figure(self._backend, FigurePlan(panels=[panel], figsize=figsize))

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
        prepared = prepare_categorical_data(
            df=df,
            category_col=category_col,
            p_col=p_col,
            category_order=category_order,
        )
        panel = categorical_spec(
            prepared,
            significance_threshold=significance_threshold,
            title=title or "Categorical Manhattan Plot",
        )
        return render_figure(self._backend, FigurePlan(panels=[panel], figsize=figsize))

    def plot_qq(
        self,
        df: pd.DataFrame,
        *,
        config: GenomeWideConfig = GenomeWideConfig(),
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
            config: Column names; only ``config.p_col`` is read.
            show_confidence_band: If True, show 95% confidence band.
            show_lambda: If True, show genomic inflation factor in title.
            figsize: Figure size as (width, height).
            title: Plot title. If None and show_lambda is True, shows lambda.

        Returns:
            Figure object (type depends on backend).

        Example:
            >>> fig = plotter.plot_qq(gwas_df, config=GenomeWideConfig(p_col="pvalue"))
        """
        qq = prepare_qq_data(df, p_col=config.p_col)
        panel = QQPanelSpec(
            qq_df=qq.frame,
            show_confidence_band=show_confidence_band,
            title=title
            or qq_title(qq.lambda_gc, show_lambda=show_lambda, compact=False),
            title_fontsize=14,
        )
        return render_figure(self._backend, FigurePlan(panels=[panel], figsize=figsize))

    def plot_manhattan_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        *,
        config: GenomeWideConfig = GenomeWideConfig(),
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
            config: Column names and chromosome order, shared by every frame.
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

        prepared = prepare_genomewide_frames(gwas_dfs, config, species=self.species)
        return render_figure(
            self._backend,
            FigurePlan(
                panels=stacked_manhattan_specs(
                    prepared,
                    significance_threshold=significance_threshold,
                    panel_labels=panel_labels,
                ),
                figsize=figsize,
                height_ratios=[figsize[1] / n_gwas] * n_gwas,
                first_panel_title=title,
                hspace=0.1,
            ),
        )

    def plot_manhattan_qq(
        self,
        df: pd.DataFrame,
        *,
        config: GenomeWideConfig = GenomeWideConfig(),
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
            config: Column names and chromosome order.
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

        manhattan = prepare_genomewide_frames([df], config, species=self.species)[0]
        qq = prepare_qq_data(df, p_col=config.p_col)
        return render_figure(
            self._backend,
            FigurePlan(
                panels=[
                    manhattan_spec(
                        manhattan,
                        significance_threshold=significance_threshold,
                        x_label="Chromosome",
                        title="Manhattan Plot",
                        title_fontsize=12,
                    ),
                    QQPanelSpec(
                        qq_df=qq.frame,
                        show_confidence_band=show_confidence_band,
                        title=qq_title(
                            qq.lambda_gc, show_lambda=show_lambda, compact=False
                        ),
                        title_fontsize=12,
                    ),
                ],
                figsize=figsize,
                n_cols=2,
                width_ratios=[2.5, 1],
                suptitle=title,
                top=0.90 if title else 0.95,
            ),
        )

    def plot_manhattan_qq_stacked(
        self,
        gwas_dfs: List[pd.DataFrame],
        *,
        config: GenomeWideConfig = GenomeWideConfig(),
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
            config: Column names and chromosome order, shared by every frame.
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

        manhattans = prepare_genomewide_frames(gwas_dfs, config, species=self.species)
        specs = stacked_manhattan_specs(
            manhattans,
            significance_threshold=significance_threshold,
            panel_labels=panel_labels,
        )
        panels = []
        for index, (spec, df) in enumerate(zip(specs, gwas_dfs)):
            qq = prepare_qq_data(df, p_col=config.p_col)
            panels.append(spec)
            panels.append(
                QQPanelSpec(
                    qq_df=qq.frame,
                    show_confidence_band=show_confidence_band,
                    title=qq_title(qq.lambda_gc, show_lambda=show_lambda, compact=True),
                    title_fontsize=10,
                    label_fontsize=10,
                    x_label="Expected $-\\log_{10}(p)$"
                    if index == n_gwas - 1
                    else None,
                )
            )
        return render_figure(
            self._backend,
            FigurePlan(
                panels=panels,
                figsize=figsize,
                n_cols=2,
                width_ratios=[2.5, 1],
                suptitle=title,
                top=0.90 if title else 0.95,
                hspace=0.15,
            ),
        )

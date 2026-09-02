"""Backend-neutral rendering modules for plot families.

The public plotters prepare data and describe the figure they want.  The
renderer owns panel composition, axes, labels, legends, and layout.  It uses
the existing primitive backend contract internally so custom backends keep
working while the rendering seam gains depth.
"""

from typing import Any, List, Optional, Sequence, Tuple

import pandas as pd

from ._manhattan_panel import (
    ManhattanPanelSpec,
    chromosome_ticks,
    manhattan_spec,
    render_manhattan_panel,
    shared_manhattan_limits,
)
from ._plotter_utils import (
    MANHATTAN_CATEGORICAL_POINT_SIZE,
    POINT_EDGE_COLOR,
    QQ_CI_ALPHA,
    QQ_CI_COLOR,
    QQ_EDGE_WIDTH,
    QQ_POINT_COLOR,
    QQ_POINT_SIZE,
    SIGNIFICANCE_LINE_COLOR,
)
from .backends.base import PlotBackend


class ManhattanQQRenderer:
    """Render prepared Manhattan and QQ figure intent through one seam.

    This is the first semantic renderer over the legacy primitive backend
    contract.  Plotters provide prepared DataFrames and presentation intent;
    this module owns the repeated figure and panel policy.
    """

    def __init__(self, backend: PlotBackend):
        self._backend = backend

    def render_manhattan(
        self,
        prepared_df: pd.DataFrame,
        *,
        figsize: Tuple[float, float],
        significance_threshold: Optional[float],
        title: Optional[str],
    ) -> Any:
        """Render one prepared Manhattan panel."""
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        render_manhattan_panel(
            self._backend,
            axes[0],
            manhattan_spec(
                prepared_df,
                significance_threshold=significance_threshold,
                x_label="Chromosome",
                title=title or "Manhattan Plot",
            ),
        )
        self._backend.finalize_layout(fig)
        return fig

    def render_categorical(
        self,
        prepared_df: pd.DataFrame,
        *,
        figsize: Tuple[float, float],
        significance_threshold: Optional[float],
        title: Optional[str],
    ) -> Any:
        """Render one prepared categorical Manhattan panel."""
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        cat_order = prepared_df.attrs["category_order"]
        category_centers = prepared_df.attrs["category_centers"]
        render_manhattan_panel(
            self._backend,
            axes[0],
            ManhattanPanelSpec(
                prepared_df=prepared_df,
                x_col="_x_pos",
                group_col="_cat_str",
                group_order=cat_order,
                x_limits=(-0.5, len(cat_order) - 0.5),
                tick_positions=[category_centers[category] for category in cat_order],
                tick_labels=cat_order,
                significance_threshold=significance_threshold,
                point_size=MANHATTAN_CATEGORICAL_POINT_SIZE,
                tick_fontsize=10,
                tick_rotation=45,
                tick_ha="right",
                x_label="Category",
                title=title or "Categorical Manhattan Plot",
            ),
        )
        self._backend.finalize_layout(fig)
        return fig

    def render_qq(
        self,
        qq_df: pd.DataFrame,
        *,
        figsize: Tuple[float, float],
        show_confidence_band: bool,
        show_lambda: bool,
        title: Optional[str],
    ) -> Any:
        """Render one prepared QQ panel."""
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
        )
        ax = axes[0]
        self._render_qq_panel(ax, qq_df, show_confidence_band)
        self._set_qq_labels_and_title(
            ax,
            qq_df,
            show_lambda=show_lambda,
            title=title,
            fontsize=12,
            title_fontsize=14,
        )
        self._backend.finalize_layout(fig)
        return fig

    def render_manhattan_stacked(
        self,
        prepared_dfs: Sequence[pd.DataFrame],
        *,
        figsize: Tuple[float, float],
        significance_threshold: Optional[float],
        panel_labels: Optional[List[str]],
        title: Optional[str],
    ) -> Any:
        """Render prepared Manhattan panels with shared x coordinates."""
        chrom_order = prepared_dfs[0].attrs["chrom_order"]
        chrom_centers = prepared_dfs[0].attrs["chrom_centers"]
        x_limits = shared_manhattan_limits(prepared_dfs)
        n_panels = len(prepared_dfs)
        fig, axes = self._backend.create_figure(
            n_panels=n_panels,
            height_ratios=[figsize[1] / n_panels] * n_panels,
            figsize=figsize,
            sharex=True,
        )
        ticks = chromosome_ticks(chrom_order, chrom_centers)

        for index, (ax, prepared_df) in enumerate(zip(axes, prepared_dfs)):
            render_manhattan_panel(
                self._backend,
                ax,
                manhattan_spec(
                    prepared_df,
                    x_limits=x_limits,
                    ticks=ticks,
                    significance_threshold=significance_threshold,
                    y_label_fontsize=10,
                    x_label="Chromosome" if index == n_panels - 1 else None,
                    panel_label=panel_labels[index]
                    if panel_labels and index < len(panel_labels)
                    else None,
                ),
            )

        if title:
            self._backend.set_title(axes[0], title, fontsize=14)
        self._backend.finalize_layout(fig, hspace=0.1)
        return fig

    def render_manhattan_qq(
        self,
        manhattan_df: pd.DataFrame,
        qq_df: pd.DataFrame,
        *,
        figsize: Tuple[float, float],
        significance_threshold: Optional[float],
        show_confidence_band: bool,
        show_lambda: bool,
        title: Optional[str],
    ) -> Any:
        """Render one side-by-side Manhattan and QQ figure."""
        fig, axes = self._backend.create_figure_grid(
            n_rows=1,
            n_cols=2,
            width_ratios=[2.5, 1],
            figsize=figsize,
        )
        render_manhattan_panel(
            self._backend,
            axes[0],
            manhattan_spec(
                manhattan_df,
                significance_threshold=significance_threshold,
                x_label="Chromosome",
                title="Manhattan Plot",
                title_fontsize=12,
            ),
        )
        self._render_qq_panel(axes[1], qq_df, show_confidence_band)
        self._set_qq_labels_and_title(
            axes[1],
            qq_df,
            show_lambda=show_lambda,
            title=None,
            fontsize=12,
            title_fontsize=12,
        )

        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90)
        else:
            self._backend.finalize_layout(fig)
        return fig

    def render_manhattan_qq_stacked(
        self,
        manhattan_dfs: Sequence[pd.DataFrame],
        qq_dfs: Sequence[pd.DataFrame],
        *,
        figsize: Tuple[float, float],
        significance_threshold: Optional[float],
        show_confidence_band: bool,
        show_lambda: bool,
        panel_labels: Optional[List[str]],
        title: Optional[str],
    ) -> Any:
        """Render stacked side-by-side Manhattan and QQ panels."""
        chrom_order = manhattan_dfs[0].attrs["chrom_order"]
        chrom_centers = manhattan_dfs[0].attrs["chrom_centers"]
        x_limits = shared_manhattan_limits(manhattan_dfs)
        n_panels = len(manhattan_dfs)
        fig, axes = self._backend.create_figure_grid(
            n_rows=n_panels,
            n_cols=2,
            width_ratios=[2.5, 1],
            figsize=figsize,
        )
        ticks = chromosome_ticks(chrom_order, chrom_centers)

        for index, (manhattan_df, qq_df) in enumerate(zip(manhattan_dfs, qq_dfs)):
            manhattan_ax = axes[index * 2]
            qq_ax = axes[index * 2 + 1]
            render_manhattan_panel(
                self._backend,
                manhattan_ax,
                manhattan_spec(
                    manhattan_df,
                    x_limits=x_limits,
                    ticks=ticks,
                    significance_threshold=significance_threshold,
                    y_label_fontsize=10,
                    x_label="Chromosome" if index == n_panels - 1 else None,
                    panel_label=panel_labels[index]
                    if panel_labels and index < len(panel_labels)
                    else None,
                ),
            )
            self._render_qq_panel(qq_ax, qq_df, show_confidence_band)
            self._set_qq_labels_and_title(
                qq_ax,
                qq_df,
                show_lambda=show_lambda,
                title=None,
                fontsize=10,
                x_label="Expected $-\\log_{10}(p)$" if index == n_panels - 1 else None,
                y_label="Observed $-\\log_{10}(p)$",
                stacked=True,
                title_fontsize=10,
            )

        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90, hspace=0.15)
        else:
            self._backend.finalize_layout(fig, hspace=0.15)
        return fig

    def _render_qq_panel(
        self,
        ax: Any,
        qq_df: pd.DataFrame,
        show_confidence_band: bool,
    ) -> None:
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
        self._backend.line(
            ax,
            x=pd.Series([0, max_val]),
            y=pd.Series([0, max_val]),
            color=SIGNIFICANCE_LINE_COLOR,
            linestyle="--",
            linewidth=1,
            zorder=2,
        )
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
        self._backend.set_xlim(ax, 0, max_val * 1.05)
        self._backend.set_ylim(ax, 0, max_val * 1.05)

    def _set_qq_labels_and_title(
        self,
        ax: Any,
        qq_df: pd.DataFrame,
        *,
        show_lambda: bool,
        title: Optional[str],
        fontsize: int,
        title_fontsize: Optional[int] = None,
        x_label: Optional[str] = r"Expected $-\log_{10}(p)$",
        y_label: str = r"Observed $-\log_{10}(p)$",
        stacked: bool = False,
    ) -> None:
        if x_label is not None:
            self._backend.set_xlabel(ax, x_label, fontsize=fontsize)
        self._backend.set_ylabel(ax, y_label, fontsize=fontsize)
        if title:
            plot_title = title
        elif show_lambda:
            lambda_gc = qq_df.attrs["lambda_gc"]
            plot_title = (
                f"λ = {lambda_gc:.3f}" if stacked else f"QQ Plot (λ = {lambda_gc:.3f})"
            )
        else:
            plot_title = "QQ" if stacked else "QQ Plot"
        self._backend.set_title(ax, plot_title, fontsize=title_fontsize or fontsize + 2)

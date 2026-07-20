"""Backend-neutral rendering modules for plot families.

The public plotters prepare data and describe the figure they want.  The
renderer owns panel composition, axes, labels, legends, and layout.  It uses
the existing primitive backend contract internally so custom backends keep
working while the rendering seam gains depth.
"""

from typing import Any, List, Optional, Sequence, Tuple

import pandas as pd

from ._plotter_utils import (
    MANHATTAN_CATEGORICAL_POINT_SIZE,
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
    QQ_CI_ALPHA,
    QQ_CI_COLOR,
    QQ_EDGE_WIDTH,
    QQ_POINT_COLOR,
    QQ_POINT_SIZE,
    SIGNIFICANCE_LINE_COLOR,
    add_significance_line,
)
from .backends.base import PlotBackend


def _padded_ymax(y_max: float) -> float:
    """Return a useful upper y-limit for Manhattan panels."""
    return max(y_max * 1.1, 1.0) if pd.notna(y_max) else 1.0


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
        self._render_manhattan_panel(
            axes[0],
            prepared_df,
            prepared_df.attrs["chrom_order"],
            prepared_df.attrs["chrom_centers"],
            significance_threshold=significance_threshold,
            x_limits=None,
            y_label_fontsize=12,
            x_label="Chromosome",
            title=title or "Manhattan Plot",
            title_fontsize=14,
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
        ax = axes[0]
        cat_order = prepared_df.attrs["category_order"]
        for category in cat_order:
            category_data = prepared_df[prepared_df["_cat_str"] == category]
            if not category_data.empty:
                self._backend.scatter(
                    ax,
                    category_data["_x_pos"],
                    category_data["_neg_log_p"],
                    colors=category_data["_color"].iloc[0],
                    sizes=MANHATTAN_CATEGORICAL_POINT_SIZE,
                    marker="o",
                    edgecolor=POINT_EDGE_COLOR,
                    linewidth=MANHATTAN_EDGE_WIDTH,
                    zorder=2,
                )

        add_significance_line(self._backend, ax, significance_threshold)
        category_centers = prepared_df.attrs["category_centers"]
        positions = [category_centers[category] for category in cat_order]
        self._backend.set_xticks(
            ax, positions, cat_order, fontsize=10, rotation=45, ha="right"
        )
        self._backend.set_xlim(ax, -0.5, len(cat_order) - 0.5)
        self._backend.set_ylim(ax, 0, _padded_ymax(prepared_df["_neg_log_p"].max()))
        self._backend.set_xlabel(ax, "Category", fontsize=12)
        self._backend.set_ylabel(ax, r"$-\log_{10}(p)$", fontsize=12)
        self._backend.set_title(ax, title or "Categorical Manhattan Plot", fontsize=14)
        self._backend.hide_spines(ax, ["top", "right"])
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
        self._backend.hide_spines(ax, ["top", "right"])
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
        x_limits = self._shared_manhattan_limits(prepared_dfs)
        n_panels = len(prepared_dfs)
        fig, axes = self._backend.create_figure(
            n_panels=n_panels,
            height_ratios=[figsize[1] / n_panels] * n_panels,
            figsize=figsize,
            sharex=True,
        )
        positions, labels = self._chromosome_ticks(chrom_order, chrom_centers)

        for index, (ax, prepared_df) in enumerate(zip(axes, prepared_dfs)):
            self._render_manhattan_panel(
                ax,
                prepared_df,
                chrom_order,
                chrom_centers,
                significance_threshold=significance_threshold,
                x_limits=x_limits,
                y_label_fontsize=10,
                x_label="Chromosome" if index == n_panels - 1 else None,
                title=None,
                panel_label=panel_labels[index]
                if panel_labels and index < len(panel_labels)
                else None,
                tick_positions=positions,
                tick_labels=labels,
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
        self._render_manhattan_panel(
            axes[0],
            manhattan_df,
            manhattan_df.attrs["chrom_order"],
            manhattan_df.attrs["chrom_centers"],
            significance_threshold=significance_threshold,
            x_limits=None,
            y_label_fontsize=12,
            x_label="Chromosome",
            title="Manhattan Plot",
            title_fontsize=12,
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
        self._backend.hide_spines(axes[1], ["top", "right"])

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
        x_limits = self._shared_manhattan_limits(manhattan_dfs)
        n_panels = len(manhattan_dfs)
        fig, axes = self._backend.create_figure_grid(
            n_rows=n_panels,
            n_cols=2,
            width_ratios=[2.5, 1],
            figsize=figsize,
        )
        positions, labels = self._chromosome_ticks(chrom_order, chrom_centers)

        for index, (manhattan_df, qq_df) in enumerate(zip(manhattan_dfs, qq_dfs)):
            manhattan_ax = axes[index * 2]
            qq_ax = axes[index * 2 + 1]
            self._render_manhattan_panel(
                manhattan_ax,
                manhattan_df,
                chrom_order,
                chrom_centers,
                significance_threshold=significance_threshold,
                x_limits=x_limits,
                y_label_fontsize=10,
                x_label="Chromosome" if index == n_panels - 1 else None,
                title=None,
                panel_label=panel_labels[index]
                if panel_labels and index < len(panel_labels)
                else None,
                tick_positions=positions,
                tick_labels=labels,
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
            self._backend.hide_spines(qq_ax, ["top", "right"])

        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90, hspace=0.15)
        else:
            self._backend.finalize_layout(fig, hspace=0.15)
        return fig

    def _render_manhattan_panel(
        self,
        ax: Any,
        prepared_df: pd.DataFrame,
        chrom_order: List[str],
        chrom_centers: dict[str, float],
        *,
        significance_threshold: Optional[float],
        x_limits: Optional[Tuple[float, float]],
        y_label_fontsize: int,
        x_label: Optional[str],
        title: Optional[str],
        panel_label: Optional[str] = None,
        tick_positions: Optional[List[float]] = None,
        tick_labels: Optional[List[str]] = None,
        title_fontsize: Optional[int] = None,
    ) -> None:
        """Apply shared Manhattan panel policy to one backend axis."""
        self._render_manhattan_points(ax, prepared_df, chrom_order)
        add_significance_line(self._backend, ax, significance_threshold)

        if x_limits is None:
            x_min = prepared_df["_cumulative_pos"].min()
            x_max = prepared_df["_cumulative_pos"].max()
            x_padding = (x_max - x_min) * 0.01
            x_limits = (x_min - x_padding, x_max + x_padding)
        self._backend.set_xlim(ax, *x_limits)
        self._backend.set_ylim(ax, 0, _padded_ymax(prepared_df["_neg_log_p"].max()))

        if tick_positions is None or tick_labels is None:
            tick_positions, tick_labels = self._chromosome_ticks(
                chrom_order, chrom_centers
            )
        self._backend.set_xticks(ax, tick_positions, tick_labels, fontsize=8)
        if x_label:
            self._backend.set_xlabel(ax, x_label, fontsize=12)
        self._backend.set_ylabel(ax, r"$-\log_{10}(p)$", fontsize=y_label_fontsize)
        if title:
            self._backend.set_title(
                ax,
                title,
                fontsize=title_fontsize or (14 if y_label_fontsize == 12 else 12),
            )
        if panel_label:
            self._backend.add_panel_label(ax, panel_label)
        self._backend.hide_spines(ax, ["top", "right"])

    def _render_manhattan_points(
        self,
        ax: Any,
        prepared_df: pd.DataFrame,
        chrom_order: List[str],
    ) -> None:
        for chrom in chrom_order:
            chrom_data = prepared_df[prepared_df["_chrom_str"] == chrom]
            if not chrom_data.empty:
                self._backend.scatter(
                    ax,
                    chrom_data["_cumulative_pos"],
                    chrom_data["_neg_log_p"],
                    colors=chrom_data["_color"].iloc[0],
                    sizes=MANHATTAN_POINT_SIZE,
                    marker="o",
                    edgecolor=POINT_EDGE_COLOR,
                    linewidth=MANHATTAN_EDGE_WIDTH,
                    zorder=2,
                )

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

    @staticmethod
    def _chromosome_ticks(
        chrom_order: List[str], chrom_centers: dict[str, float]
    ) -> Tuple[List[float], List[str]]:
        labels = [chrom for chrom in chrom_order if chrom in chrom_centers]
        return [chrom_centers[chrom] for chrom in labels], labels

    @staticmethod
    def _shared_manhattan_limits(
        prepared_dfs: Sequence[pd.DataFrame],
    ) -> Tuple[float, float]:
        x_min = min(df["_cumulative_pos"].min() for df in prepared_dfs)
        x_max = max(df["_cumulative_pos"].max() for df in prepared_dfs)
        x_padding = (x_max - x_min) * 0.01
        return x_min - x_padding, x_max + x_padding

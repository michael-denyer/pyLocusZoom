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
    manhattan_spec,
    render_manhattan_panel,
)
from ._plotter_utils import MANHATTAN_CATEGORICAL_POINT_SIZE
from ._qq_panel import QQPanelSpec, qq_title, render_qq_panel
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
            height_ratios=[1.0],
            figsize=figsize,
        )
        render_manhattan_panel(
            self._backend,
            axes[0],
            ManhattanPanelSpec(
                prepared_df=prepared_df,
                x_col="_x_pos",
                group_col="_cat_str",
                layout=prepared_df.attrs["layout"],
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
            height_ratios=[1.0],
            figsize=figsize,
        )
        render_qq_panel(
            self._backend,
            axes[0],
            QQPanelSpec(
                qq_df=qq_df,
                show_confidence_band=show_confidence_band,
                title=title
                or qq_title(
                    qq_df.attrs["lambda_gc"], show_lambda=show_lambda, compact=False
                ),
                title_fontsize=14,
            ),
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
        n_panels = len(prepared_dfs)
        fig, axes = self._backend.create_figure(
            height_ratios=[figsize[1] / n_panels] * n_panels,
            figsize=figsize,
            sharex=True,
        )
        for ax, spec in zip(
            axes,
            self._stacked_manhattan_specs(
                prepared_dfs,
                significance_threshold=significance_threshold,
                panel_labels=panel_labels,
            ),
        ):
            render_manhattan_panel(self._backend, ax, spec)

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
        render_qq_panel(
            self._backend,
            axes[1],
            QQPanelSpec(
                qq_df=qq_df,
                show_confidence_band=show_confidence_band,
                title=qq_title(
                    qq_df.attrs["lambda_gc"], show_lambda=show_lambda, compact=False
                ),
                title_fontsize=12,
            ),
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
        n_panels = len(manhattan_dfs)
        fig, axes = self._backend.create_figure_grid(
            n_rows=n_panels,
            n_cols=2,
            width_ratios=[2.5, 1],
            figsize=figsize,
        )
        specs = self._stacked_manhattan_specs(
            manhattan_dfs,
            significance_threshold=significance_threshold,
            panel_labels=panel_labels,
        )

        for index, (spec, qq_df) in enumerate(zip(specs, qq_dfs)):
            render_manhattan_panel(self._backend, axes[index * 2], spec)
            render_qq_panel(
                self._backend,
                axes[index * 2 + 1],
                QQPanelSpec(
                    qq_df=qq_df,
                    show_confidence_band=show_confidence_band,
                    title=qq_title(
                        qq_df.attrs["lambda_gc"], show_lambda=show_lambda, compact=True
                    ),
                    title_fontsize=10,
                    label_fontsize=10,
                    x_label="Expected $-\\log_{10}(p)$"
                    if index == n_panels - 1
                    else None,
                ),
            )

        if title:
            self._backend.set_suptitle(fig, title, fontsize=14)
            self._backend.finalize_layout(fig, top=0.90, hspace=0.15)
        else:
            self._backend.finalize_layout(fig, hspace=0.15)
        return fig

    @staticmethod
    def _stacked_manhattan_specs(
        prepared_dfs: Sequence[pd.DataFrame],
        *,
        significance_threshold: Optional[float],
        panel_labels: Optional[List[str]],
    ) -> List[ManhattanPanelSpec]:
        """Build specs for vertically stacked panels sharing one genome layout.

        Only the bottom panel carries the x-axis label, since the panels share
        one x axis.

        Args:
            prepared_dfs: Frames from ``prepare_manhattan_frames``, top to
                bottom.
            significance_threshold: P-value for the significance line, or None.
            panel_labels: Corner label per panel, or None.

        Returns:
            One spec per frame, in the same order.
        """
        n_panels = len(prepared_dfs)
        return [
            manhattan_spec(
                prepared_df,
                significance_threshold=significance_threshold,
                y_label_fontsize=10,
                x_label="Chromosome" if index == n_panels - 1 else None,
                panel_label=panel_labels[index]
                if panel_labels and index < len(panel_labels)
                else None,
            )
            for index, prepared_df in enumerate(prepared_dfs)
        ]

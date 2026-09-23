"""One Manhattan panel: the typed request, which draws itself.

A Manhattan panel, a categorical (PheWAS-style) panel, and each half of a
Miami plot are the same nine drawing steps over different data columns, tick
sets, and limits. :class:`ManhattanPanelSpec` names those differences so the
Manhattan plotter's plans and the Miami plan share one policy instead of
three copies.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, TypeVar

import numpy as np
import pandas as pd

from ..backends.base import PlotBackend
from ..backends.hover import HoverConfig, HoverDataBuilder
from ..config import GenomeWideStyle
from ..exceptions import ValidationError
from ..manhattan import PreparedManhattan
from ._shared import (
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
    SUGGESTIVE_LINE_COLOR,
    add_significance_line,
)

T = TypeVar("T")


def styled(override: Optional[T], default: T) -> T:
    """Return a style field the caller set, or the panel's own default."""
    return default if override is None else override


def padded_ymax(y_max: float, headroom: float) -> float:
    """Return a useful upper y-limit for a Manhattan panel."""
    return max(y_max * (1 + headroom), 1.0) if pd.notna(y_max) else 1.0


@dataclass(frozen=True)
class ManhattanPanelSpec:
    """One Manhattan-style panel's data and presentation policy.

    Attributes:
        prepared: The frame from ``prepare_genomewide_frames`` or
            ``prepare_categorical_data``, with the x and group columns its
            preparation created and the layout every panel of the figure
            shares.
        significance_threshold: P-value to draw the significance line at, or
            None to draw no line.
        suggestive_threshold: P-value to draw the suggestive line at, or
            None to draw no line.
        point_size: Marker size.
        tick_fontsize: X tick label size.
        tick_rotation: X tick label rotation.
        tick_ha: X tick label horizontal alignment.
        x_label: X axis label, or None for none.
        y_label_fontsize: Y axis label size.
        title: Panel title, or None for none.
        title_fontsize: Panel title size.
        panel_label: Corner label, or None for none.
        panel_label_y_frac: Fractional height of the corner label.
        invert_y: Draw the y axis descending, as the lower Miami panel does.
        hover: Hover column mapping, or None for no tooltips.
        style: Caller styling. A field it sets overrides the matching
            field above.
    """

    prepared: PreparedManhattan
    significance_threshold: Optional[float] = None
    suggestive_threshold: Optional[float] = None
    point_size: int = MANHATTAN_POINT_SIZE
    tick_fontsize: int = 8
    tick_rotation: int = 0
    tick_ha: str = "center"
    x_label: Optional[str] = None
    y_label_fontsize: int = 12
    title: Optional[str] = None
    title_fontsize: int = 14
    panel_label: Optional[str] = None
    panel_label_y_frac: float = 0.95
    invert_y: bool = False
    hover: Optional[HoverConfig] = None
    style: GenomeWideStyle = GenomeWideStyle()

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw this panel onto a backend axis."""
        df, layout = self.prepared.frame, self.prepared.layout
        style = self.style
        for group in layout.order:
            group_data = df[df[self.prepared.group_col] == group]
            if group_data.empty:
                continue
            hover_data = (
                HoverDataBuilder(self.hover).build(group_data)
                if self.hover is not None
                else None
            )
            backend.scatter(
                ax,
                group_data[self.prepared.x_col],
                group_data["neglog10p"],
                colors=group_data["_color"].iloc[0],
                sizes=styled(style.point_size, self.point_size),
                marker="o",
                edgecolor=POINT_EDGE_COLOR,
                linewidth=styled(style.point_edge_width, MANHATTAN_EDGE_WIDTH),
                zorder=2,
                hover_data=hover_data,
                alpha=style.point_alpha,
            )

        line_kwargs = dict(linestyle=style.line_style, linewidth=style.line_width)
        add_significance_line(backend, ax, self.significance_threshold, **line_kwargs)
        add_significance_line(
            backend,
            ax,
            self.suggestive_threshold,
            color=SUGGESTIVE_LINE_COLOR,
            **line_kwargs,
        )
        backend.set_xlim(ax, *layout.x_limits)
        line_levels = [
            -np.log10(threshold)
            for threshold in (self.significance_threshold, self.suggestive_threshold)
            if threshold is not None
        ]
        y_max = padded_ymax(
            max([df["neglog10p"].max(), *line_levels]), style.y_headroom
        )
        if self.invert_y:
            backend.set_ylim(ax, y_max, 0)
        else:
            backend.set_ylim(ax, 0, y_max)
        backend.set_xticks(
            ax,
            layout.tick_positions[:: style.tick_step],
            layout.tick_labels[:: style.tick_step],
            fontsize=styled(style.tick_label_fontsize, self.tick_fontsize),
            rotation=styled(style.tick_rotation, self.tick_rotation),
            ha=self.tick_ha,
        )
        if style.tick_label_fontsize is not None:
            backend.set_tick_fontsize(ax, style.tick_label_fontsize)
        if self.x_label:
            backend.set_xlabel(
                ax, self.x_label, fontsize=styled(style.axis_label_fontsize, 12)
            )
        backend.set_ylabel(
            ax,
            r"$-\log_{10}(p)$",
            fontsize=styled(style.axis_label_fontsize, self.y_label_fontsize),
        )
        if self.title:
            backend.set_title(
                ax,
                self.title,
                fontsize=styled(style.panel_title_fontsize, self.title_fontsize),
                fontweight=style.title_fontweight,
            )
        if self.panel_label:
            backend.add_panel_label(
                ax, self.panel_label, y_frac=self.panel_label_y_frac
            )


def stacked_manhattan_specs(
    prepared: Sequence[PreparedManhattan],
    *,
    significance_threshold: Optional[float],
    panel_labels: Optional[Sequence[str]],
    style: GenomeWideStyle = GenomeWideStyle(),
) -> List[ManhattanPanelSpec]:
    """Build specs for vertically stacked panels sharing one genome layout.

    Only the bottom panel carries the x-axis label, since the panels share
    one x axis.

    Args:
        prepared: Values from ``prepare_manhattan_frames``, top to bottom.
        significance_threshold: P-value for the significance line, or None.
        panel_labels: Corner label per panel, or None.
        style: Caller styling, shared by every panel.

    Returns:
        One spec per frame, in the same order.

    Raises:
        ValidationError: If ``panel_labels`` does not hold one label per frame.
    """
    n_panels = len(prepared)
    if panel_labels is not None and len(panel_labels) != n_panels:
        raise ValidationError(
            f"panel_labels length ({len(panel_labels)}) must match "
            f"number of GWAS DataFrames ({n_panels})"
        )
    return [
        ManhattanPanelSpec(
            value,
            significance_threshold=significance_threshold,
            y_label_fontsize=10,
            x_label="Chromosome" if index == n_panels - 1 else None,
            panel_label=panel_labels[index] if panel_labels is not None else None,
            style=style,
        )
        for index, value in enumerate(prepared)
    ]

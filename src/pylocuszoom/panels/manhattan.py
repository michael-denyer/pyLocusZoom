"""One Manhattan panel: the typed request and the function that draws it.

A Manhattan panel, a categorical (PheWAS-style) panel, and each half of a
Miami plot are the same nine drawing steps over different data columns, tick
sets, and limits. :class:`ManhattanPanelSpec` names those differences so the
Manhattan plotter's plans and the Miami plan share one policy instead of
three copies.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence

import pandas as pd

from .._plotter_utils import (
    MANHATTAN_CATEGORICAL_POINT_SIZE,
    MANHATTAN_EDGE_WIDTH,
    MANHATTAN_POINT_SIZE,
    POINT_EDGE_COLOR,
    add_significance_line,
)
from ..backends.base import PlotBackend
from ..backends.hover import HoverConfig, HoverDataBuilder
from ..manhattan import PanelLayout, PreparedManhattan


def padded_ymax(y_max: float) -> float:
    """Return a useful upper y-limit for a Manhattan panel."""
    return max(y_max * 1.1, 1.0) if pd.notna(y_max) else 1.0


@dataclass(frozen=True)
class ManhattanPanelSpec:
    """One Manhattan-style panel's data and presentation policy.

    Attributes:
        prepared_df: The ``frame`` of a :class:`~.manhattan.PreparedManhattan`.
        x_col: Column holding each point's x coordinate.
        group_col: Column the scatter loop groups by, one colour per group.
        layout: Group order, x limits, and ticks, shared by every panel of
            the figure.
        significance_threshold: P-value to draw the significance line at, or
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
        hover: Hover column mapping, or None for no tooltips. Built only
            when the backend reports ``supports_hover``, since matplotlib
            discards the frame the builder would allocate per group.
    """

    prepared_df: pd.DataFrame
    x_col: str
    group_col: str
    layout: PanelLayout
    significance_threshold: Optional[float] = None
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

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw this panel onto a backend axis."""
        render_manhattan_panel(backend, ax, self)


def manhattan_spec(
    prepared: PreparedManhattan,
    *,
    significance_threshold: Optional[float] = None,
    x_label: Optional[str] = None,
    y_label_fontsize: int = 12,
    title: Optional[str] = None,
    title_fontsize: int = 14,
    panel_label: Optional[str] = None,
    panel_label_y_frac: float = 0.95,
    invert_y: bool = False,
    hover: Optional[HoverConfig] = None,
) -> ManhattanPanelSpec:
    """Build a genomic-position panel spec from a prepared Manhattan frame.

    The keyword arguments are the :class:`ManhattanPanelSpec` fields the
    single, stacked and Miami call sites vary; the tick styling fields keep
    the spec's own defaults. ``test_manhattan_spec_defaults_match_the_spec``
    fails if the two lists of defaults drift.

    Args:
        prepared: One value from ``prepare_manhattan_frames``, carrying the
            frame and the shared :class:`~.manhattan.GenomeLayout`.
        significance_threshold: P-value to draw the significance line at, or
            None to draw no line.
        x_label: X axis label, or None for none.
        y_label_fontsize: Y axis label size.
        title: Panel title, or None for none.
        title_fontsize: Panel title size.
        panel_label: Corner label, or None for none.
        panel_label_y_frac: Fractional height of the corner label.
        invert_y: Draw the y axis descending, as the lower Miami panel does.
        hover: Hover column mapping, or None for no tooltips.

    Returns:
        The panel spec.
    """
    return ManhattanPanelSpec(
        prepared_df=prepared.frame,
        x_col="_cumulative_pos",
        group_col="_chrom_str",
        layout=prepared.layout,
        significance_threshold=significance_threshold,
        x_label=x_label,
        y_label_fontsize=y_label_fontsize,
        title=title,
        title_fontsize=title_fontsize,
        panel_label=panel_label,
        panel_label_y_frac=panel_label_y_frac,
        invert_y=invert_y,
        hover=hover,
    )


def categorical_spec(
    prepared: PreparedManhattan,
    *,
    significance_threshold: Optional[float],
    title: str,
) -> ManhattanPanelSpec:
    """Build a category-axis panel spec from a prepared categorical frame.

    Args:
        prepared: The value from ``prepare_categorical_data``.
        significance_threshold: P-value to draw the significance line at, or
            None to draw no line.
        title: Panel title.

    Returns:
        The panel spec, with the larger points and rotated ticks a category
        axis needs.
    """
    return ManhattanPanelSpec(
        prepared_df=prepared.frame,
        x_col="_x_pos",
        group_col="_cat_str",
        layout=prepared.layout,
        significance_threshold=significance_threshold,
        point_size=MANHATTAN_CATEGORICAL_POINT_SIZE,
        tick_fontsize=10,
        tick_rotation=45,
        tick_ha="right",
        x_label="Category",
        title=title,
    )


def stacked_manhattan_specs(
    prepared: Sequence[PreparedManhattan],
    *,
    significance_threshold: Optional[float],
    panel_labels: Optional[Sequence[str]],
) -> List[ManhattanPanelSpec]:
    """Build specs for vertically stacked panels sharing one genome layout.

    Only the bottom panel carries the x-axis label, since the panels share
    one x axis.

    Args:
        prepared: Values from ``prepare_manhattan_frames``, top to bottom.
        significance_threshold: P-value for the significance line, or None.
        panel_labels: Corner label per panel, or None.

    Returns:
        One spec per frame, in the same order.
    """
    n_panels = len(prepared)
    return [
        manhattan_spec(
            value,
            significance_threshold=significance_threshold,
            y_label_fontsize=10,
            x_label="Chromosome" if index == n_panels - 1 else None,
            panel_label=panel_labels[index]
            if panel_labels and index < len(panel_labels)
            else None,
        )
        for index, value in enumerate(prepared)
    ]


def render_manhattan_panel(
    backend: PlotBackend, ax: Any, spec: ManhattanPanelSpec
) -> None:
    """Draw one Manhattan-style panel onto a backend axis.

    Args:
        backend: Backend that owns the drawing primitives.
        ax: Axis to draw onto.
        spec: The panel's data and presentation policy.
    """
    df = spec.prepared_df
    for group in spec.layout.order:
        group_data = df[df[spec.group_col] == group]
        if group_data.empty:
            continue
        hover_data = None
        if spec.hover is not None and backend.supports_hover:
            hover_data = HoverDataBuilder(spec.hover).build_dataframe(group_data)
        backend.scatter(
            ax,
            group_data[spec.x_col],
            group_data["neglog10p"],
            colors=group_data["_color"].iloc[0],
            sizes=spec.point_size,
            marker="o",
            edgecolor=POINT_EDGE_COLOR,
            linewidth=MANHATTAN_EDGE_WIDTH,
            zorder=2,
            hover_data=hover_data,
        )

    add_significance_line(backend, ax, spec.significance_threshold)
    backend.set_xlim(ax, *spec.layout.x_limits)
    y_max = padded_ymax(df["neglog10p"].max())
    if spec.invert_y:
        backend.set_ylim(ax, y_max, 0)
    else:
        backend.set_ylim(ax, 0, y_max)
    backend.set_xticks(
        ax,
        spec.layout.tick_positions,
        spec.layout.tick_labels,
        fontsize=spec.tick_fontsize,
        rotation=spec.tick_rotation,
        ha=spec.tick_ha,
    )
    if spec.x_label:
        backend.set_xlabel(ax, spec.x_label, fontsize=12)
    backend.set_ylabel(ax, r"$-\log_{10}(p)$", fontsize=spec.y_label_fontsize)
    if spec.title:
        backend.set_title(ax, spec.title, fontsize=spec.title_fontsize)
    if spec.panel_label:
        backend.add_panel_label(ax, spec.panel_label, y_frac=spec.panel_label_y_frac)

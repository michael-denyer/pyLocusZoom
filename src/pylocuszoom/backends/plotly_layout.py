"""Plotly subplot geometry: the panel value type and pure layout helpers.

Plotly addresses a subplot by a linear index and positions a legend in paper
coordinates against an axis domain. Both are arithmetic over a figure rather
than drawing, so they live beside the backend as plain functions and one value
type. A reader following a protocol method resolves nothing through the
instance to reach them.
"""

from typing import Any, List, NamedTuple, Optional, Tuple

import plotly.graph_objects as go

# Matplotlib legend `loc` vocabulary mapped to Plotly (xanchor, yanchor).
# "best" has no Plotly equivalent, so it takes the upper-right default.
_LEGEND_ANCHORS = {
    "best": ("right", "top"),
    "upper right": ("right", "top"),
    "upper left": ("left", "top"),
    "upper center": ("center", "top"),
    "lower right": ("right", "bottom"),
    "lower left": ("left", "bottom"),
    "lower center": ("center", "bottom"),
    "center right": ("right", "middle"),
    "center left": ("left", "middle"),
    "right": ("right", "middle"),
    "center": ("center", "middle"),
}
_LEGEND_X = {"left": 0.01, "center": 0.5, "right": 0.99}


class _Panel(NamedTuple):
    """One subplot of a Plotly figure, resolved from a renderer's panel handle.

    Plotly names axes by a linear subplot index: subplot (1,1) is ``xaxis`` and
    ``x``, subplot (1,2) is ``xaxis2`` and ``x2``. Owning that arithmetic here
    keeps the naming rule in one place instead of at every call site.
    """

    fig: go.Figure
    row: int
    col: int = 1
    n_cols: int = 1

    @property
    def subplot_idx(self) -> int:
        """Plotly's linear index for this subplot, 1-based."""
        return (self.row - 1) * self.n_cols + self.col

    def axis(self, kind: str) -> str:
        """Layout key for this subplot's axis, such as ``"xaxis3"``.

        Args:
            kind: Either ``"xaxis"`` or ``"yaxis"``.

        Returns:
            The layout key, unsuffixed for the first subplot.
        """
        idx = self.subplot_idx
        return f"{kind}{idx}" if idx > 1 else kind

    def secondary_ref(self) -> str:
        """Trace-level reference for this subplot's secondary y-axis.

        Offset by 100 so the name cannot collide with the primary axes, which
        Plotly numbers yaxis, yaxis2, ..., yaxisN for N subplots. That supports
        up to 99 subplot rows.

        Returns:
            The secondary axis reference, such as ``"y100"``.
        """
        return f"y{100 + self.subplot_idx - 1}"

    def ref(self, kind: str) -> str:
        """Trace-level reference for this subplot's axis, such as ``"x3"``.

        Args:
            kind: Either ``"x"`` or ``"y"``.

        Returns:
            The axis reference, unsuffixed for the first subplot.
        """
        idx = self.subplot_idx
        return f"{kind}{idx}" if idx > 1 else kind

    @property
    def xref(self) -> str:
        """Trace-level reference for this subplot's x-axis."""
        return self.ref("x")

    @property
    def yref(self) -> str:
        """Trace-level reference for this subplot's y-axis."""
        return self.ref("y")


class _SecondaryAxis(NamedTuple):
    """A panel's secondary y-axis, as returned by ``create_twin_axis``.

    Carries the panel it overlays, so a drawing primitive reaches the figure
    and the panel's x-axis through the handle it is given.
    """

    panel: _Panel
    yref: str

    @property
    def fig(self) -> go.Figure:
        """The figure the panel belongs to."""
        return self.panel.fig

    @property
    def xref(self) -> str:
        """Trace-level reference for the panel's x-axis, which is shared."""
        return self.panel.xref


def secondary_axis_key(secondary_ref: str) -> str:
    """Layout key for a secondary axis given its trace-level reference.

    Args:
        secondary_ref: A reference such as ``"y100"``.

    Returns:
        The matching layout key, such as ``"yaxis100"``.
    """
    if secondary_ref.startswith("y"):
        return "yaxis" + secondary_ref[1:]
    return secondary_ref


def panel_y(fig: go.Figure, row: int, vertical: str) -> float:
    """Find a y-coordinate in paper coords inside a subplot row's domain.

    Args:
        fig: The figure carrying the subplot.
        row: Row number, 1-indexed.
        vertical: One of ``"bottom"``, ``"middle"``, or ``"top"``.

    Returns:
        The y-coordinate in paper coordinates.
    """
    yaxis = getattr(fig.layout, _Panel(fig, row).axis("yaxis"), None)
    domain = yaxis.domain if yaxis and yaxis.domain else (0.01, 0.99)
    if vertical == "bottom":
        return domain[0]
    if vertical == "middle":
        return (domain[0] + domain[1]) / 2
    return domain[1]


def configure_legend(
    fig: go.Figure, row: int, legend_key: str, title: str, loc: str
) -> None:
    """Position and style one of a figure's legends against a panel row.

    Args:
        fig: The figure carrying the legend.
        legend_key: Layout key for this legend, such as ``"legend2"``.
        row: Row number the legend belongs to, 1-indexed.
        title: Legend title, already in display form.
        loc: Matplotlib legend location vocabulary, such as ``"upper right"``.
    """
    horizontal, vertical = _LEGEND_ANCHORS.get(loc, _LEGEND_ANCHORS["upper right"])
    fig.update_layout(
        **{
            legend_key: dict(
                title=dict(text=title),
                x=_LEGEND_X[horizontal],
                y=panel_y(fig, row, vertical),
                xanchor=horizontal,
                yanchor=vertical,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        }
    )


def x_range(panel: _Panel, xaxis_name: str) -> Optional[Tuple[float, float]]:
    """Resolve a panel's x-range, falling back to the extent of trace data.

    Args:
        panel: The panel whose range is wanted.
        xaxis_name: Layout key for the panel's x-axis.

    Returns:
        The (min, max) range, or None when the figure carries no x data.
    """
    xaxis = getattr(panel.fig.layout, xaxis_name, None)
    if xaxis and xaxis.range:
        return xaxis.range
    x_vals: List[Any] = []
    for trace in panel.fig.data:
        if getattr(trace, "x", None) is not None:
            x_vals.extend([v for v in trace.x if v is not None])
    return (min(x_vals), max(x_vals)) if x_vals else None

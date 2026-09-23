"""Plotly subplot geometry: the panel value type and pure layout helpers.

Plotly addresses a subplot by a linear index and positions a legend in paper
coordinates against an axis domain. Both are arithmetic over a figure rather
than drawing, so they live beside the backend as plain functions and one value
type. A reader following a protocol method resolves nothing through the
instance to reach them.
"""

from typing import NamedTuple, Optional, Tuple

import plotly.graph_objects as go


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


def panel_top(panel: _Panel) -> float:
    """The top of a subplot's domain, in paper coordinates.

    Args:
        panel: The subplot.

    Returns:
        The y-coordinate in paper coordinates.
    """
    yaxis = getattr(panel.fig.layout, panel.axis("yaxis"), None)
    return yaxis.domain[1] if yaxis and yaxis.domain else 0.99


def configure_legend(panel: _Panel, legend_key: str, title: str) -> None:
    """Anchor one of a figure's legends in a panel's upper-right corner.

    Args:
        panel: The subplot the legend belongs to.
        legend_key: Layout key for this legend, such as ``"legend2"``.
        title: Legend title, already in display form.
    """
    panel.fig.update_layout(
        **{
            legend_key: dict(
                title=dict(text=title),
                x=0.99,
                y=panel_top(panel),
                xanchor="right",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
            )
        }
    )


def x_range(panel: _Panel, xaxis_name: str) -> Optional[Tuple[float, float]]:
    """A panel's x-axis limits, as set with ``set_xlim``.

    A shared-x figure links its panels with ``matches``, so every axis in
    the group spans whatever any one of them was given. The traces are not
    consulted: a panel whose data is narrower than the region must still be
    labelled for the region.

    Args:
        panel: The panel whose range is wanted.
        xaxis_name: Layout key for the panel's x-axis.

    Returns:
        The (min, max) range, or None when no axis in the group has limits.
    """
    layout = panel.fig.layout

    def root(name: str) -> str:
        seen = set()
        while name not in seen:
            seen.add(name)
            axis = getattr(layout, name, None)
            if axis is None or not axis.matches:
                return name
            name = axis.matches.replace("x", "xaxis", 1)
        return name

    group = [xaxis_name] + [
        name
        for name in layout
        if name.startswith("xaxis")
        and name != xaxis_name
        and root(name) == root(xaxis_name)
    ]
    for name in group:
        axis = getattr(layout, name, None)
        if axis is not None and axis.range:
            return axis.range
    return None

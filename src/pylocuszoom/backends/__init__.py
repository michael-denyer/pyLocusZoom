"""Pluggable plotting backends for pyLocusZoom.

Supports matplotlib (default), plotly, and bokeh backends.
"""

from typing import Literal

from .base import PlotBackend
from .matplotlib_backend import MatplotlibBackend

BackendType = Literal["matplotlib", "plotly", "bokeh"]

_BACKENDS: dict[str, type[PlotBackend]] = {
    "matplotlib": MatplotlibBackend,
}


def get_backend(name: BackendType) -> PlotBackend:
    """Get a backend instance by name.

    Args:
        name: Backend name ('matplotlib', 'plotly', or 'bokeh').

    Returns:
        Instantiated backend.

    Raises:
        ValueError: If backend name is invalid.
        ImportError: If backend dependencies are not installed.
    """
    if name == "plotly":
        from .plotly_backend import PlotlyBackend

        _BACKENDS["plotly"] = PlotlyBackend
    elif name == "bokeh":
        from .bokeh_backend import BokehBackend

        _BACKENDS["bokeh"] = BokehBackend

    if name not in _BACKENDS:
        raise ValueError(
            f"Unknown backend: {name}. Available: matplotlib, plotly, bokeh"
        )

    return _BACKENDS[name]()


__all__ = ["PlotBackend", "BackendType", "get_backend", "MatplotlibBackend"]

"""Pluggable plotting backends for pyLocusZoom.

Supports matplotlib (default), plotly, and bokeh backends.
"""

from typing import Literal

from .base import PlotBackend
from .matplotlib_backend import MatplotlibBackend

BackendType = Literal["matplotlib", "plotly", "bokeh"]

# LaTeX to Unicode conversion table for interactive backends
_LATEX_TO_UNICODE = [
    (r"$-\log_{10}$ P", "-log₁₀(P)"),
    (r"$-\log_{10}$", "-log₁₀"),
    (r"\log_{10}", "log₁₀"),
    (r"$r^2$", "r²"),
    (r"$R^2$", "R²"),
]


def convert_latex_to_unicode(label: str) -> str:
    """Convert LaTeX-style labels to Unicode for display in interactive backends.

    Args:
        label: Label text possibly containing LaTeX notation.

    Returns:
        Label with LaTeX converted to Unicode characters.
    """
    for latex, unicode_str in _LATEX_TO_UNICODE:
        if latex in label:
            label = label.replace(latex, unicode_str)
    return label.replace("$", "")


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


__all__ = [
    "PlotBackend",
    "BackendType",
    "get_backend",
    "MatplotlibBackend",
    "convert_latex_to_unicode",
]

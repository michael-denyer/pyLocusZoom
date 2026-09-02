"""Pluggable plotting backends for pyLocusZoom.

Supports matplotlib (default), plotly, and bokeh backends.

Backend Registration:
    Backends are registered using the @register_backend decorator.
    This enables extensibility without modifying core code.

    @register_backend("custom")
    class CustomBackend:
        ...

Error Behavior:
    When an optional backend (plotly, bokeh) is requested but not installed,
    get_backend() raises ImportError with installation instructions.
"""

from typing import Literal, get_args

from .base import Mappable, PlotBackend, SupportsSNPLabels

BackendType = Literal["matplotlib", "plotly", "bokeh"]

BUILTIN_BACKENDS: tuple[str, ...] = get_args(BackendType)
"""Every backend name shipped with pyLocusZoom, in registration order."""

# LaTeX to Unicode conversion table for interactive backends
_LATEX_TO_UNICODE = [
    (r"$-\log_{10}$ P", "-log10(P)"),
    (r"$-\log_{10}$", "-log10"),
    (r"\log_{10}", "log10"),
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


# Backend registry - populated by @register_backend decorator
_BACKENDS: dict[str, type[PlotBackend]] = {}


def register_backend(name: str):
    """Decorator to register a backend class.

    Args:
        name: Backend name (e.g., "matplotlib", "plotly", "bokeh").

    Returns:
        Decorator function that registers the class.

    Example:
        @register_backend("custom")
        class CustomBackend:
            ...
    """

    def decorator(cls: type[PlotBackend]) -> type[PlotBackend]:
        _BACKENDS[name] = cls
        return cls

    return decorator


def get_backend(name: BackendType) -> PlotBackend:
    """Get a backend instance by name.

    Args:
        name: Backend name ('matplotlib', 'plotly', or 'bokeh').

    Returns:
        Instantiated backend.

    Raises:
        ImportError: If plotly or bokeh is not installed.
        ValueError: If backend name is completely unknown.
    """
    # Ensure matplotlib is always registered (it's always available)
    if "matplotlib" not in _BACKENDS:
        from .matplotlib_backend import MatplotlibBackend  # noqa: F401

    # Try lazy import for optional backends
    if name not in _BACKENDS:
        if name == "plotly":
            try:
                from .plotly_backend import PlotlyBackend  # noqa: F401
            except ImportError:
                raise ImportError(
                    "Plotly is required for the plotly backend. "
                    "Install it with: pip install plotly"
                )
        elif name == "bokeh":
            try:
                from .bokeh_backend import BokehBackend  # noqa: F401
            except ImportError:
                raise ImportError(
                    "Bokeh is required for the bokeh backend. "
                    "Install it with: pip install bokeh"
                )

    if name not in _BACKENDS:
        available = list(_BACKENDS.keys())
        raise ValueError(f"Unknown backend: {name}. Available: {available}")

    return _BACKENDS[name]()


# Lazy import for backward compatibility - MatplotlibBackend available at module level
# The actual registration happens when matplotlib_backend is imported
def __getattr__(name: str):
    """Lazy attribute access for backward compatibility."""
    if name == "MatplotlibBackend":
        from .matplotlib_backend import MatplotlibBackend

        return MatplotlibBackend
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "Mappable",
    "PlotBackend",
    "SupportsSNPLabels",
    "BackendType",
    "BUILTIN_BACKENDS",
    "get_backend",
    "register_backend",
    "MatplotlibBackend",
    "convert_latex_to_unicode",
]

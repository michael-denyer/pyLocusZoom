"""Logging switches for pylocuszoom.

pyLocusZoom logs through loguru and is silent until asked, following loguru's
pattern for libraries: importing the package calls
``logger.disable("pylocuszoom")`` and adds or removes no handler, so the host
application's sinks are exactly as it left them. ``enable_logging`` is the one
switch that turns pyLocusZoom's records on.

Log records are diagnostics. A layer the plot has to leave out (the gene
track, the recombination overlay or LD colouring) is reported as a
``UserWarning`` instead, so it reaches the user whether logging is on or not.

Usage:
    >>> from pylocuszoom.logging import enable_logging, disable_logging
    >>> enable_logging("DEBUG")  # DEBUG level for troubleshooting
    >>> disable_logging()  # silent again
"""

import sys
from typing import Optional

from loguru import logger

logger.disable("pylocuszoom")

_FORMAT = "<level>{level: <8}</level> | <cyan>pylocuszoom</cyan> | {message}"
_handler_id: Optional[int] = None


def _remove_handler() -> None:
    """Remove the sink enable_logging added, if another caller has not already."""
    global _handler_id
    if _handler_id is not None:
        try:
            logger.remove(_handler_id)
        except ValueError:
            pass  # A global logger.remove() elsewhere already took it.
        _handler_id = None


def enable_logging(level: str = "INFO", sink=sys.stderr) -> None:
    """Send pyLocusZoom's records at ``level`` and above to ``sink``.

    A second call replaces the sink the first one added. The records also
    reach any sink the application added without a filter, loguru's default
    stderr sink included; an application that routes everything through its
    own sinks can call loguru's ``logger.enable("pylocuszoom")`` instead.

    Args:
        level: Log level ("DEBUG", "INFO", "WARNING", "ERROR").
        sink: Output destination (default: stderr).

    Example:
        >>> from pylocuszoom.logging import enable_logging
        >>> enable_logging()  # INFO level
        >>> enable_logging("DEBUG")  # DEBUG level for troubleshooting
    """
    global _handler_id
    _remove_handler()
    _handler_id = logger.add(sink, level=level, format=_FORMAT, filter="pylocuszoom")
    logger.enable("pylocuszoom")


def disable_logging() -> None:
    """Silence pyLocusZoom's records and remove the sink enable_logging added."""
    _remove_handler()
    logger.disable("pylocuszoom")

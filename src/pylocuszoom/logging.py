"""Logging configuration for pylocuszoom.

Provides logging with sensible defaults:
- Logging is enabled by default at INFO level
- Uses loguru (a hard dependency)
- Users can adjust level via enable_logging() or disable via disable_logging()

Usage:
    >>> from pylocuszoom.logging import enable_logging, disable_logging
    >>> enable_logging("DEBUG")  # Enable DEBUG level for troubleshooting
    >>> disable_logging()  # Suppress all logging output
"""

import sys
from typing import Optional

from loguru import logger as _loguru_logger

_LOGURU_DEFAULT_HANDLER_ID = 0


class _LoguruWrapper:
    """Wrapper around loguru logger with enable/disable support.

    Errors always emit (even when logging is disabled) so that callers
    using logger.error() for real failures are never silently swallowed.
    """

    def __init__(self):
        self._enabled = False
        self._handler_id = None
        self._error_handler_id = None
        self._remove_handler(_LOGURU_DEFAULT_HANDLER_ID)

    def _add_error_handler(self) -> None:
        """Add error-only handler (used when main handler is disabled)."""
        if self._error_handler_id is None:
            self._error_handler_id = _loguru_logger.add(
                sys.stderr,
                level="ERROR",
                format="<level>{level: <8}</level> | <cyan>pylocuszoom</cyan> | {message}",
                filter=lambda record: record["name"].startswith("pylocuszoom"),
            )

    def _remove_handler(self, handler_id: Optional[int]) -> None:
        """Remove a handler, tolerating an id another module already removed.

        A global ``loguru.remove()`` elsewhere invalidates the ids stored here,
        so a stale id is an expected state rather than a fault.
        """
        if handler_id is None:
            return
        try:
            _loguru_logger.remove(handler_id)
        except ValueError:
            pass

    def enable(self, level: str = "INFO", sink=sys.stderr) -> None:
        """Enable logging at the specified level."""
        # Remove error-only handler (main handler covers errors)
        self._remove_handler(self._error_handler_id)
        self._error_handler_id = None
        self._remove_handler(self._handler_id)
        self._handler_id = _loguru_logger.add(
            sink,
            level=level,
            format="<level>{level: <8}</level> | <cyan>pylocuszoom</cyan> | {message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )
        self._enabled = True

    def disable(self) -> None:
        """Disable logging (errors still emit via error-only handler)."""
        self._remove_handler(self._handler_id)
        self._handler_id = None
        # Re-add error-only handler so errors still reach stderr
        self._add_error_handler()
        self._enabled = False

    def debug(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            _loguru_logger.opt(depth=1).debug(msg, *args, **kwargs)

    def info(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            _loguru_logger.opt(depth=1).info(msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            _loguru_logger.opt(depth=1).warning(msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs) -> None:
        _loguru_logger.opt(depth=1).error(msg, *args, **kwargs)


# Create the logger instance
logger = _LoguruWrapper()

# Enable logging at INFO level by default
logger.enable("INFO")


def enable_logging(level: str = "INFO", sink=sys.stderr) -> None:
    """Enable logging output.

    Args:
        level: Log level ("DEBUG", "INFO", "WARNING", "ERROR").
        sink: Output destination (default: stderr).

    Example:
        >>> from pylocuszoom.logging import enable_logging
        >>> enable_logging()  # INFO level
        >>> enable_logging("DEBUG")  # DEBUG level for troubleshooting
    """
    logger.enable(level, sink)


def disable_logging() -> None:
    """Disable logging output."""
    logger.disable()

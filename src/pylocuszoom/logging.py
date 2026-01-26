"""Logging configuration for pylocuszoom.

Provides logging with sensible defaults:
- Logging is enabled by default at INFO level
- Uses loguru (included as dependency)
- Users can adjust level via enable_logging() or disable via disable_logging()

Usage:
    >>> from pylocuszoom.logging import enable_logging, disable_logging
    >>> enable_logging("DEBUG")  # Enable DEBUG level for troubleshooting
    >>> disable_logging()  # Suppress all logging output
"""

import sys

# Try to use loguru, fall back to stdlib logging
try:
    from loguru import logger as _loguru_logger

    _HAS_LOGURU = True
except ImportError:
    import logging as _stdlib_logging

    _HAS_LOGURU = False


class _LoguruWrapper:
    """Wrapper around loguru logger with enable/disable support."""

    def __init__(self):
        self._enabled = False
        self._handler_id = None
        # Remove default handler
        _loguru_logger.remove()

    def enable(self, level: str = "INFO", sink=sys.stderr) -> None:
        """Enable logging at the specified level."""
        if self._handler_id is not None:
            try:
                _loguru_logger.remove(self._handler_id)
            except ValueError:
                # Handler was already removed (e.g., by another module calling logger.remove())
                pass
        self._handler_id = _loguru_logger.add(
            sink,
            level=level,
            format="<level>{level: <8}</level> | <cyan>pylocuszoom</cyan> | {message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )
        self._enabled = True

    def disable(self) -> None:
        """Disable logging."""
        if self._handler_id is not None:
            try:
                _loguru_logger.remove(self._handler_id)
            except ValueError:
                # Handler was already removed (e.g., by another module calling logger.remove())
                pass
            self._handler_id = None
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
        if self._enabled:
            _loguru_logger.opt(depth=1).error(msg, *args, **kwargs)


class _StdlibWrapper:
    """Wrapper around stdlib logging with enable/disable support."""

    def __init__(self):
        self._logger = _stdlib_logging.getLogger("pylocuszoom")
        self._logger.setLevel(_stdlib_logging.WARNING)
        self._handler = None
        self._enabled = False

    def enable(self, level: str = "INFO", sink=sys.stderr) -> None:
        """Enable logging at the specified level."""
        if self._handler is not None:
            self._logger.removeHandler(self._handler)
        self._handler = _stdlib_logging.StreamHandler(sink)
        self._handler.setFormatter(
            _stdlib_logging.Formatter("%(levelname)-8s | pylocuszoom | %(message)s")
        )
        self._logger.addHandler(self._handler)
        self._logger.setLevel(getattr(_stdlib_logging, level.upper()))
        self._enabled = True

    def disable(self) -> None:
        """Disable logging."""
        if self._handler is not None:
            self._logger.removeHandler(self._handler)
            self._handler = None
        self._logger.setLevel(_stdlib_logging.WARNING)
        self._enabled = False

    def debug(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            self._logger.debug(msg, *args, **kwargs)

    def info(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            self._logger.info(msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            self._logger.warning(msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs) -> None:
        if self._enabled:
            self._logger.error(msg, *args, **kwargs)


# Create the logger instance
if _HAS_LOGURU:
    logger = _LoguruWrapper()
else:
    logger = _StdlibWrapper()

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

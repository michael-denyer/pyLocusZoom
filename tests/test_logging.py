"""Tests for logging utilities."""

from loguru import logger as _loguru_logger


class TestLoggingWrapper:
    """Tests for the logging wrapper."""

    def test_enable_after_external_handler_removal(self):
        """enable_logging should not raise when handler was removed externally.

        This can happen when another module (e.g., utils/__init__.py) calls
        logger.remove() globally, invalidating handler IDs stored by the wrapper.
        """
        from snp_scope_plot.logging import enable_logging

        # Simulate another module removing all handlers
        _loguru_logger.remove()

        # Should not raise ValueError
        enable_logging("INFO")

    def test_disable_after_external_handler_removal(self):
        """disable_logging should not raise when handler was removed externally."""
        from snp_scope_plot.logging import disable_logging, enable_logging

        # First enable to get a handler ID stored
        enable_logging("INFO")

        # Simulate another module removing all handlers
        _loguru_logger.remove()

        # Should not raise ValueError
        disable_logging()

    def test_enable_disable_cycle(self):
        """Enable and disable should work in sequence without errors."""
        from snp_scope_plot.logging import disable_logging, enable_logging

        enable_logging("DEBUG")
        disable_logging()
        enable_logging("INFO")
        disable_logging()

    def test_multiple_enables_without_disable(self):
        """Multiple enable calls should not accumulate handlers."""
        from snp_scope_plot.logging import disable_logging, enable_logging

        enable_logging("DEBUG")
        enable_logging("INFO")
        enable_logging("WARNING")
        disable_logging()

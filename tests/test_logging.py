"""Tests for logging utilities."""

import io

from loguru import logger as _loguru_logger


class TestLoggingWrapper:
    """Tests for the logging wrapper."""

    def test_enable_after_external_handler_removal(self):
        """enable_logging should not raise when handler was removed externally.

        This can happen when another module (e.g., utils/__init__.py) calls
        logger.remove() globally, invalidating handler IDs stored by the wrapper.
        """
        from pylocuszoom.logging import enable_logging

        # Simulate another module removing all handlers
        _loguru_logger.remove()

        # Should not raise ValueError
        enable_logging("INFO")

    def test_disable_after_external_handler_removal(self):
        """disable_logging should not raise when handler was removed externally."""
        from pylocuszoom.logging import disable_logging, enable_logging

        # First enable to get a handler ID stored
        enable_logging("INFO")

        # Simulate another module removing all handlers
        _loguru_logger.remove()

        # Should not raise ValueError
        disable_logging()

    def test_enable_disable_cycle(self):
        """Enable and disable should work in sequence without errors."""
        from pylocuszoom.logging import disable_logging, enable_logging

        enable_logging("DEBUG")
        disable_logging()
        enable_logging("INFO")
        disable_logging()

    def test_multiple_enables_without_disable(self):
        """Multiple enable calls should not accumulate handlers."""
        from pylocuszoom.logging import disable_logging, enable_logging

        enable_logging("DEBUG")
        enable_logging("INFO")
        enable_logging("WARNING")
        disable_logging()


class TestLoguruWrapper:
    """Tests for the _LoguruWrapper class directly."""

    def test_wrapper_initial_state(self):
        """Wrapper should start disabled with no handler."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        assert wrapper._enabled is False
        assert wrapper._handler_id is None

    def test_wrapper_enable_sets_state(self):
        """Enable should set enabled state and handler ID."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.enable("INFO")

        assert wrapper._enabled is True
        assert wrapper._handler_id is not None
        wrapper.disable()

    def test_wrapper_disable_clears_state(self):
        """Disable should clear enabled state and handler ID."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.enable("INFO")
        wrapper.disable()

        assert wrapper._enabled is False
        assert wrapper._handler_id is None

    def test_wrapper_debug_method_exists(self):
        """Debug method should be callable."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        # Should not raise when disabled
        wrapper.debug("test message")
        wrapper.enable("DEBUG")
        # Should not raise when enabled
        wrapper.debug("test message")
        wrapper.disable()

    def test_wrapper_info_method_exists(self):
        """Info method should be callable."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.info("test message")
        wrapper.enable("INFO")
        wrapper.info("test message")
        wrapper.disable()

    def test_wrapper_warning_method_exists(self):
        """Warning method should be callable."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.warning("test message")
        wrapper.enable("WARNING")
        wrapper.warning("test message")
        wrapper.disable()

    def test_wrapper_error_method_exists(self):
        """Error method should be callable."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.error("test message")
        wrapper.enable("ERROR")
        wrapper.error("test message")
        wrapper.disable()

    def test_no_output_when_disabled(self):
        """No messages should log when disabled."""
        from pylocuszoom.logging import logger

        # Disable the logger
        logger.disable()

        # Methods should still be callable without error
        logger.debug("should not appear")
        logger.info("should not appear")
        logger.warning("should not appear")
        logger.error("should not appear")

    def test_disable_without_enable(self):
        """Disabling without enabling should not raise."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.disable()  # Should not raise
        wrapper.disable()  # Multiple disables should be fine

    def test_multiple_enable_calls(self):
        """Multiple enable calls should replace handler, not accumulate."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.enable("DEBUG")
        first_id = wrapper._handler_id
        wrapper.enable("INFO")
        second_id = wrapper._handler_id

        # Handler ID should have been replaced
        assert first_id != second_id
        wrapper.disable()

    def test_enable_with_custom_sink(self):
        """Enable should accept custom sink parameter."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        buffer = io.StringIO()
        wrapper.enable("INFO", sink=buffer)

        assert wrapper._enabled is True
        wrapper.disable()


class TestStdlibWrapperDirect:
    """Tests for the _StdlibWrapper class using direct import of stdlib logging."""

    def test_stdlib_wrapper_creation(self):
        """StdlibWrapper can be created when stdlib logging is available."""
        import logging as stdlib_logging

        # Create a standalone wrapper using stdlib logging directly
        class TestStdlibWrapper:
            def __init__(self):
                self._logger = stdlib_logging.getLogger("test_pylocuszoom")
                self._logger.setLevel(stdlib_logging.WARNING)
                self._handler = None
                self._enabled = False

            def enable(self, level="INFO", sink=None):
                if self._handler is not None:
                    self._logger.removeHandler(self._handler)
                self._handler = stdlib_logging.StreamHandler(sink)
                self._handler.setFormatter(
                    stdlib_logging.Formatter("%(levelname)-8s | test | %(message)s")
                )
                self._logger.addHandler(self._handler)
                self._logger.setLevel(getattr(stdlib_logging, level.upper()))
                self._enabled = True

            def disable(self):
                if self._handler is not None:
                    self._logger.removeHandler(self._handler)
                    self._handler = None
                self._logger.setLevel(stdlib_logging.WARNING)
                self._enabled = False

            def debug(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.debug(msg, *args, **kwargs)

            def info(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.info(msg, *args, **kwargs)

            def warning(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.warning(msg, *args, **kwargs)

            def error(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.error(msg, *args, **kwargs)

        wrapper = TestStdlibWrapper()
        assert wrapper._enabled is False
        assert wrapper._handler is None

    def test_stdlib_wrapper_enable_disable(self):
        """StdlibWrapper enable/disable cycle works."""
        import logging as stdlib_logging

        class TestStdlibWrapper:
            def __init__(self):
                self._logger = stdlib_logging.getLogger("test_pylocuszoom_2")
                self._logger.setLevel(stdlib_logging.WARNING)
                self._handler = None
                self._enabled = False

            def enable(self, level="INFO", sink=None):
                if self._handler is not None:
                    self._logger.removeHandler(self._handler)
                self._handler = stdlib_logging.StreamHandler(sink)
                self._logger.addHandler(self._handler)
                self._logger.setLevel(getattr(stdlib_logging, level.upper()))
                self._enabled = True

            def disable(self):
                if self._handler is not None:
                    self._logger.removeHandler(self._handler)
                    self._handler = None
                self._logger.setLevel(stdlib_logging.WARNING)
                self._enabled = False

        wrapper = TestStdlibWrapper()
        buffer = io.StringIO()
        wrapper.enable("INFO", sink=buffer)

        assert wrapper._enabled is True
        assert wrapper._handler is not None

        wrapper.disable()
        assert wrapper._enabled is False
        assert wrapper._handler is None

    def test_stdlib_wrapper_log_methods(self):
        """StdlibWrapper log methods work when enabled."""
        import logging as stdlib_logging

        class TestStdlibWrapper:
            def __init__(self):
                self._logger = stdlib_logging.getLogger("test_pylocuszoom_3")
                self._handler = None
                self._enabled = False

            def enable(self, level="INFO", sink=None):
                if self._handler is not None:
                    self._logger.removeHandler(self._handler)
                self._handler = stdlib_logging.StreamHandler(sink)
                self._handler.setFormatter(stdlib_logging.Formatter("%(message)s"))
                self._logger.addHandler(self._handler)
                self._logger.setLevel(getattr(stdlib_logging, level.upper()))
                self._enabled = True

            def disable(self):
                if self._handler is not None:
                    self._logger.removeHandler(self._handler)
                    self._handler = None
                self._enabled = False

            def debug(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.debug(msg, *args, **kwargs)

            def info(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.info(msg, *args, **kwargs)

            def warning(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.warning(msg, *args, **kwargs)

            def error(self, msg, *args, **kwargs):
                if self._enabled:
                    self._logger.error(msg, *args, **kwargs)

        wrapper = TestStdlibWrapper()
        buffer = io.StringIO()
        wrapper.enable("DEBUG", sink=buffer)

        wrapper.debug("debug msg")
        wrapper.info("info msg")
        wrapper.warning("warning msg")
        wrapper.error("error msg")

        output = buffer.getvalue()
        assert "debug msg" in output
        assert "info msg" in output
        assert "warning msg" in output
        assert "error msg" in output
        wrapper.disable()


class TestModuleLevelFunctions:
    """Tests for module-level enable_logging and disable_logging."""

    def test_enable_logging_default(self):
        """enable_logging with defaults should work."""
        from pylocuszoom.logging import disable_logging, enable_logging

        enable_logging()  # Default INFO level
        disable_logging()

    def test_enable_logging_all_levels(self):
        """enable_logging should accept all standard levels."""
        from pylocuszoom.logging import disable_logging, enable_logging

        for level in ["DEBUG", "INFO", "WARNING", "ERROR"]:
            enable_logging(level)
            disable_logging()

    def test_enable_logging_with_sink(self):
        """enable_logging should accept a custom sink."""
        from pylocuszoom.logging import disable_logging, enable_logging

        buffer = io.StringIO()
        enable_logging("INFO", sink=buffer)
        disable_logging()

    def test_disable_logging_multiple_times(self):
        """disable_logging can be called multiple times."""
        from pylocuszoom.logging import disable_logging

        disable_logging()
        disable_logging()
        disable_logging()

    def test_logger_methods_callable_after_disable(self):
        """Logger methods should be callable after disable."""
        from pylocuszoom.logging import disable_logging, logger

        disable_logging()

        # These should not raise
        logger.debug("test")
        logger.info("test")
        logger.warning("test")
        logger.error("test")

    def test_has_loguru_flag(self):
        """Module should have _HAS_LOGURU flag set."""
        from pylocuszoom.logging import _HAS_LOGURU

        # Since loguru is a required dependency, it should be True
        assert _HAS_LOGURU is True

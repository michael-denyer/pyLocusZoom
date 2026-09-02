"""Tests for logging utilities and exception hierarchy."""

import io
import subprocess
import sys

import pandas as pd
import pytest
from loguru import logger as _loguru_logger

from pylocuszoom.exceptions import (
    DataDownloadError,
    EQTLValidationError,
    FinemappingValidationError,
    ForestValidationError,
    LoaderValidationError,
    PheWASValidationError,
    PlinkError,
    PyLocusZoomError,
    ValidationError,
)


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

    def test_logger_importable_from_package(self):
        """Logger should be importable from the package top level."""
        from pylocuszoom import disable_logging, enable_logging

        enable_logging("WARNING")
        disable_logging()

    def test_logger_uses_loguru(self):
        """Logger instance should be a loguru wrapper (not stdlib)."""
        from pylocuszoom.logging import _LoguruWrapper, logger

        assert isinstance(logger, _LoguruWrapper)

    def test_error_always_emits_when_disabled(self):
        """Error messages should always be emitted even when logging disabled."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        # Disabled by default - error should still be callable without raising
        assert wrapper._enabled is False
        wrapper.error("critical failure")  # Must not raise


# =============================================================================
# Exception hierarchy tests
# =============================================================================


class TestExceptionHierarchy:
    """Tests for the pyLocusZoom exception class hierarchy."""

    def test_all_exceptions_inherit_from_pylocuszoom_error(self):
        """Every custom exception inherits from PyLocusZoomError."""
        exception_classes = [
            ValidationError,
            EQTLValidationError,
            FinemappingValidationError,
            LoaderValidationError,
            PheWASValidationError,
            ForestValidationError,
            PlinkError,
            DataDownloadError,
        ]
        for cls in exception_classes:
            assert issubclass(cls, PyLocusZoomError), (
                f"{cls.__name__} should inherit from PyLocusZoomError"
            )

    def test_validation_errors_catchable_as_valueerror(self):
        """All ValidationError subclasses can be caught with except ValueError."""
        validation_classes = [
            ValidationError,
            EQTLValidationError,
            FinemappingValidationError,
            LoaderValidationError,
            PheWASValidationError,
            ForestValidationError,
        ]
        for cls in validation_classes:
            with pytest.raises(ValueError):
                raise cls("test error")

    def test_plink_error_catchable_as_runtime_error(self):
        """PlinkError can be caught with except RuntimeError."""
        with pytest.raises(RuntimeError):
            raise PlinkError("plink failed")

    def test_data_download_error_catchable_as_runtime_error(self):
        """DataDownloadError can be caught with except RuntimeError."""
        with pytest.raises(RuntimeError):
            raise DataDownloadError("download failed")

    def test_phewas_validation_error_is_validation_error(self):
        """PheWASValidationError is a subclass of ValidationError."""
        assert issubclass(PheWASValidationError, ValidationError)

    def test_forest_validation_error_is_validation_error(self):
        """ForestValidationError is a subclass of ValidationError."""
        assert issubclass(ForestValidationError, ValidationError)

    def test_catch_all_with_pylocuszoom_error(self):
        """All library exceptions are catchable with except PyLocusZoomError."""
        for cls in [
            ValidationError,
            PheWASValidationError,
            ForestValidationError,
            PlinkError,
            DataDownloadError,
        ]:
            with pytest.raises(PyLocusZoomError):
                raise cls("test")

    def test_exception_message_preserved(self):
        """Exception message is accessible via str()."""
        msg = "Column 'pos' is missing"
        err = LoaderValidationError(msg)
        assert str(err) == msg

    def test_exceptions_importable_from_package(self):
        """All exceptions are importable from the top-level package."""
        import pylocuszoom

        assert hasattr(pylocuszoom, "PyLocusZoomError")
        assert hasattr(pylocuszoom, "PheWASValidationError")
        assert hasattr(pylocuszoom, "ForestValidationError")
        assert hasattr(pylocuszoom, "PlinkError")
        assert hasattr(pylocuszoom, "DataDownloadError")


class TestSpecializedExceptionsInUse:
    """Tests that phewas.py and forest.py raise specialized exceptions."""

    def test_phewas_validation_raises_phewas_error(self):
        """validate_phewas_df raises PheWASValidationError, not generic ValidationError."""
        from pylocuszoom.schemas import validate_phewas_df

        df = pd.DataFrame({"wrong_col": [1]})
        with pytest.raises(PheWASValidationError):
            validate_phewas_df(df)

    def test_phewas_error_also_caught_as_validation_error(self):
        """PheWAS error is still catchable as ValidationError (backward compat)."""
        from pylocuszoom.schemas import validate_phewas_df

        df = pd.DataFrame({"wrong_col": [1]})
        with pytest.raises(ValidationError):
            validate_phewas_df(df)

    def test_forest_validation_raises_forest_error(self):
        """validate_forest_df raises ForestValidationError, not generic ValidationError."""
        from pylocuszoom.schemas import validate_forest_df

        df = pd.DataFrame({"wrong_col": [1]})
        with pytest.raises(ForestValidationError):
            validate_forest_df(df)

    def test_forest_error_also_caught_as_validation_error(self):
        """Forest error is still catchable as ValidationError (backward compat)."""
        from pylocuszoom.schemas import validate_forest_df

        df = pd.DataFrame({"wrong_col": [1]})
        with pytest.raises(ValidationError):
            validate_forest_df(df)


class TestImportSideEffects:
    """Importing the package must leave the host application's sinks alone."""

    def test_host_sink_survives_import(self):
        script = (
            "import io, sys\n"
            "from loguru import logger\n"
            "buf = io.StringIO()\n"
            "logger.add(buf, format='{message}')\n"
            "import pylocuszoom.logging\n"
            "logger.info('host sink alive')\n"
            "sys.stdout.write(buf.getvalue())\n"
        )
        result = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, check=True
        )

        assert "host sink alive" in result.stdout

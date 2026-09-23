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


@pytest.fixture(autouse=True)
def restore_logging_state():
    """Put the global logger back the way this file found it.

    Most tests here end in ``disable_logging()``, which mutates process-wide
    state. Under ``pytest-randomly`` that leaked into any later test asserting
    on a warning, which then passed or failed by execution order.
    """
    from pylocuszoom.logging import disable_logging, enable_logging, logger

    was_enabled = logger._enabled
    yield
    if was_enabled:
        enable_logging("INFO")
    else:
        disable_logging()


def _emit(level, message):
    """Log ``message`` the way library code does, from a pylocuszoom module.

    The library's handlers only take records whose module is ``pylocuszoom*``,
    so a call made from this test module would be filtered out whatever the
    wrapper did.
    """
    from pylocuszoom.logging import logger

    exec(
        f"logger.{level}(message)",
        {"__name__": "pylocuszoom.probe", "logger": logger, "message": message},
    )


class TestLoggingWrapper:
    """Tests for the logging wrapper."""

    def test_enable_after_external_handler_removal(self):
        """enable_logging still attaches a working sink after a global remove().

        This can happen when another module (e.g., utils/__init__.py) calls
        logger.remove() globally, invalidating handler IDs stored by the wrapper.
        """
        from pylocuszoom.logging import enable_logging

        _loguru_logger.remove()
        sink = io.StringIO()
        enable_logging("INFO", sink=sink)

        _emit("info", "after removal")
        assert "after removal" in sink.getvalue()

    def test_disable_after_external_handler_removal(self):
        """disable_logging tolerates a handler another module already removed."""
        from pylocuszoom.logging import disable_logging, enable_logging

        sink = io.StringIO()
        enable_logging("INFO", sink=sink)
        _loguru_logger.remove()

        disable_logging()

        _emit("info", "while disabled")
        assert sink.getvalue() == ""

    def test_enable_disable_cycle(self):
        """Each enable takes its own level; each disable silences the last sink."""
        from pylocuszoom.logging import disable_logging, enable_logging

        first, second = io.StringIO(), io.StringIO()
        enable_logging("DEBUG", sink=first)
        _emit("debug", "debug one")
        disable_logging()
        _emit("info", "while disabled")
        enable_logging("INFO", sink=second)
        _emit("debug", "debug two")
        _emit("info", "info two")
        disable_logging()

        assert "debug one" in first.getvalue()
        assert "while disabled" not in first.getvalue() + second.getvalue()
        assert "debug two" not in second.getvalue()
        assert "info two" in second.getvalue()

    def test_multiple_enables_without_disable(self):
        """Multiple enable calls should not accumulate handlers."""
        from pylocuszoom.logging import enable_logging

        sinks = [io.StringIO() for _ in range(3)]
        for level, sink in zip(["DEBUG", "INFO", "WARNING"], sinks):
            enable_logging(level, sink=sink)

        _emit("warning", "only the last sink")

        assert [sink.getvalue().count("only the last sink") for sink in sinks] == [
            0,
            0,
            1,
        ]


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

    @pytest.mark.parametrize("level", ["debug", "info", "warning", "error"])
    def test_each_level_reaches_the_sink_only_while_enabled(self, level):
        """Every level method writes to the sink when enabled, and only then."""
        from pylocuszoom.logging import disable_logging, enable_logging

        sink = io.StringIO()
        enable_logging("DEBUG", sink=sink)
        _emit(level, "while enabled")
        disable_logging()
        _emit(level, "while disabled")

        assert "while enabled" in sink.getvalue()
        assert "while disabled" not in sink.getvalue()

    def test_disable_without_enable(self):
        """Disabling a wrapper that was never enabled leaves it disabled."""
        from pylocuszoom.logging import _LoguruWrapper

        wrapper = _LoguruWrapper()
        wrapper.disable()
        wrapper.disable()

        assert wrapper._enabled is False
        assert wrapper._handler_id is None

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
        """Enable routes library messages to the sink it is given."""
        from pylocuszoom.logging import enable_logging

        buffer = io.StringIO()
        enable_logging("INFO", sink=buffer)
        _emit("info", "to the custom sink")

        assert "to the custom sink" in buffer.getvalue()


class TestModuleLevelFunctions:
    """Tests for module-level enable_logging and disable_logging."""

    def test_enable_logging_defaults_to_info(self):
        """enable_logging() without a level passes INFO and above, not DEBUG."""
        from pylocuszoom.logging import enable_logging

        sink = io.StringIO()
        enable_logging(sink=sink)
        _emit("debug", "default debug")
        _emit("info", "default info")

        assert "default info" in sink.getvalue()
        assert "default debug" not in sink.getvalue()

    @pytest.mark.parametrize(
        ("level", "shown"),
        [
            ("DEBUG", ["debug", "info", "warning", "error"]),
            ("INFO", ["info", "warning", "error"]),
            ("WARNING", ["warning", "error"]),
            ("ERROR", ["error"]),
        ],
    )
    def test_enable_logging_level_is_the_floor(self, level, shown):
        """enable_logging(level) passes that level and above, nothing below."""
        from pylocuszoom.logging import enable_logging

        sink = io.StringIO()
        enable_logging(level, sink=sink)
        for method in ["debug", "info", "warning", "error"]:
            _emit(method, f"<{method}>")

        written = [
            m
            for m in ["debug", "info", "warning", "error"]
            if f"<{m}>" in sink.getvalue()
        ]
        assert written == shown

    def test_disable_logging_multiple_times(self):
        """Repeated disable_logging calls keep the library silent."""
        from pylocuszoom.logging import disable_logging, enable_logging

        sink = io.StringIO()
        enable_logging("DEBUG", sink=sink)
        disable_logging()
        disable_logging()
        disable_logging()
        _emit("warning", "after disables")

        assert sink.getvalue() == ""

    def test_logger_importable_from_package(self):
        """The package top level re-exports the logging switches."""
        import pylocuszoom
        from pylocuszoom import logging as plz_logging

        assert pylocuszoom.enable_logging is plz_logging.enable_logging
        assert pylocuszoom.disable_logging is plz_logging.disable_logging

    def test_logger_uses_loguru(self):
        """Logger instance should be a loguru wrapper (not stdlib)."""
        from pylocuszoom.logging import _LoguruWrapper, logger

        assert isinstance(logger, _LoguruWrapper)

    def test_error_always_emits_when_disabled(self, capsys):
        """Error messages reach stderr even while logging is disabled."""
        from pylocuszoom.logging import disable_logging

        disable_logging()
        _emit("error", "critical failure")

        assert "critical failure" in capsys.readouterr().err


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

    def test_default_sink_survives_import(self):
        """loguru's own stderr sink keeps printing the host's messages."""
        script = (
            "from loguru import logger\n"
            "import pylocuszoom\n"
            "logger.info('host default sink alive')\n"
            "exec(\n"
            "    'from pylocuszoom.logging import logger\\n'\n"
            "    'logger.info(\"library message\")',\n"
            "    {'__name__': 'pylocuszoom.probe'},\n"
            ")\n"
        )
        result = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, check=True
        )

        assert "host default sink alive" in result.stderr
        assert result.stderr.count("library message") == 1

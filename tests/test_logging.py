"""Tests for the logging switches."""

import io
import subprocess
import sys

import pytest
from loguru import logger as _loguru_logger


@pytest.fixture(autouse=True)
def restore_logging_state():
    """Leave logging off, the package default, whatever the test switched on.

    Under ``pytest-randomly`` a test that left logging on leaked into any later
    test asserting on output, which then passed or failed by execution order.
    """
    from pylocuszoom.logging import disable_logging

    yield
    disable_logging()


def _emit(level, message):
    """Log ``message`` the way library code does, from a pylocuszoom module.

    The library's records are those whose module is ``pylocuszoom*``, so a
    call made from this test module would be filtered out whatever the
    switches did.
    """
    from pylocuszoom.logging import logger

    exec(
        f"logger.{level}(message)",
        {"__name__": "pylocuszoom.probe", "logger": logger, "message": message},
    )


def _run(script):
    return subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    )


class TestSwitches:
    def test_enable_after_external_handler_removal(self):
        """enable_logging still attaches a working sink after a global remove()."""
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

    @pytest.mark.parametrize("level", ["debug", "info", "warning", "error"])
    def test_each_level_reaches_the_sink_only_while_enabled(self, level):
        """Every level writes to the sink when enabled, and none after disable."""
        from pylocuszoom.logging import disable_logging, enable_logging

        sink = io.StringIO()
        enable_logging("DEBUG", sink=sink)
        _emit(level, "while enabled")
        disable_logging()
        _emit(level, "while disabled")

        assert "while enabled" in sink.getvalue()
        assert "while disabled" not in sink.getvalue()

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

    def test_the_sink_takes_only_pylocuszoom_records(self):
        from pylocuszoom.logging import enable_logging

        sink = io.StringIO()
        enable_logging("DEBUG", sink=sink)
        _loguru_logger.info("host message")

        assert "host message" not in sink.getvalue()

    def test_switches_are_exported_from_the_package(self):
        import pylocuszoom
        from pylocuszoom import logging as plz_logging

        assert pylocuszoom.enable_logging is plz_logging.enable_logging
        assert pylocuszoom.disable_logging is plz_logging.disable_logging

    def test_library_modules_log_through_loguru(self):
        from pylocuszoom.logging import logger

        assert logger is _loguru_logger


class TestImportSideEffects:
    """Importing the package must leave the host application's logging alone."""

    def test_import_leaves_the_handlers_unchanged(self):
        result = _run(
            "from loguru import logger\n"
            "before = dict(logger._core.handlers)\n"
            "import pylocuszoom\n"
            "print(dict(logger._core.handlers) == before)\n"
        )

        assert result.stdout.strip() == "True"

    def test_logging_is_off_until_enabled(self):
        """The default stderr sink prints the host's messages and none of ours."""
        result = _run(
            "from loguru import logger\n"
            "import pylocuszoom\n"
            "logger.info('host default sink alive')\n"
            "exec(\n"
            "    'from pylocuszoom.logging import logger\\n'\n"
            "    'logger.warning(\"library message\")',\n"
            "    {'__name__': 'pylocuszoom.probe'},\n"
            ")\n"
        )

        assert "host default sink alive" in result.stderr
        assert "library message" not in result.stderr

    def test_constructing_a_plotter_does_not_switch_logging_on(self):
        """LocusZoomPlotter used to call enable_logging("INFO") by default."""
        result = _run(
            "import io\n"
            "from loguru import logger\n"
            "logger.remove()\n"
            "buf = io.StringIO()\n"
            "logger.add(buf, format='{message}')\n"
            "from pylocuszoom import LocusZoomPlotter\n"
            "LocusZoomPlotter(species=None)\n"
            "exec(\n"
            "    'from pylocuszoom.logging import logger\\n'\n"
            "    'logger.warning(\"library message\")',\n"
            "    {'__name__': 'pylocuszoom.probe'},\n"
            ")\n"
            "print(repr(buf.getvalue()))\n"
        )

        assert result.stdout.strip() == "''"

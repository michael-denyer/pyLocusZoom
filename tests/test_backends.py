"""Tests for backend registration and fallback."""

import warnings
from unittest.mock import patch

import pytest


class TestRegisterBackend:
    """Tests for the @register_backend decorator."""

    def test_register_backend_adds_to_registry(self):
        """@register_backend decorator adds class to _BACKENDS dict."""
        from pylocuszoom.backends import _BACKENDS, register_backend

        @register_backend("test_dummy")
        class DummyBackend:
            pass

        assert "test_dummy" in _BACKENDS
        assert _BACKENDS["test_dummy"] is DummyBackend

        # Clean up
        del _BACKENDS["test_dummy"]

    def test_register_backend_returns_class_unchanged(self):
        """Decorator returns the class unchanged."""
        from pylocuszoom.backends import _BACKENDS, register_backend

        @register_backend("test_unchanged")
        class OriginalBackend:
            def method(self):
                return "original"

        assert OriginalBackend().method() == "original"

        # Clean up
        del _BACKENDS["test_unchanged"]


class TestGetBackend:
    """Tests for get_backend function."""

    def test_get_backend_matplotlib_always_works(self):
        """get_backend('matplotlib') returns MatplotlibBackend instance."""
        from pylocuszoom.backends import get_backend
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = get_backend("matplotlib")
        assert isinstance(backend, MatplotlibBackend)

    def test_get_backend_returns_new_instance(self):
        """get_backend returns a new instance each call."""
        from pylocuszoom.backends import get_backend

        backend1 = get_backend("matplotlib")
        backend2 = get_backend("matplotlib")
        assert backend1 is not backend2

    def test_get_backend_unknown_raises_valueerror(self):
        """get_backend raises ValueError for unknown backend names."""
        from pylocuszoom.backends import get_backend

        with pytest.raises(ValueError) as exc_info:
            get_backend("nonexistent_backend")

        error_msg = str(exc_info.value)
        assert "Unknown backend" in error_msg
        assert "nonexistent_backend" in error_msg
        # Should list available backends
        assert "matplotlib" in error_msg

    def test_get_backend_plotly_works_when_installed(self):
        """get_backend('plotly') returns PlotlyBackend when plotly is available."""
        pytest.importorskip("plotly")
        from pylocuszoom.backends import get_backend
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = get_backend("plotly")
        assert isinstance(backend, PlotlyBackend)

    def test_get_backend_bokeh_works_when_installed(self):
        """get_backend('bokeh') returns BokehBackend when bokeh is available."""
        pytest.importorskip("bokeh")
        from pylocuszoom.backends import get_backend
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = get_backend("bokeh")
        assert isinstance(backend, BokehBackend)


class TestGracefulFallback:
    """Tests for fallback behavior when optional backends unavailable."""

    def test_plotly_fallback_to_matplotlib(self):
        """Falls back to matplotlib with warning when plotly import fails."""
        from pylocuszoom.backends import _BACKENDS
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        # Clear plotly from registry if present
        _BACKENDS.pop("plotly", None)

        # Mock the plotly_backend import to raise ImportError
        with patch.dict(
            "sys.modules",
            {"plotly": None, "plotly.graph_objects": None, "plotly.subplots": None},
        ):
            # Need to remove the cached import of plotly_backend
            import sys

            sys.modules.pop("pylocuszoom.backends.plotly_backend", None)

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                # Re-import get_backend to get fresh function
                from importlib import reload

                import pylocuszoom.backends

                reload(pylocuszoom.backends)
                from pylocuszoom.backends import get_backend

                backend = get_backend("plotly")

                # Should get matplotlib backend
                assert isinstance(backend, MatplotlibBackend)

                # Should have warned
                assert len(w) >= 1
                warning_messages = [str(warning.message) for warning in w]
                plotly_warnings = [
                    msg
                    for msg in warning_messages
                    if "plotly" in msg.lower() or "Plotly" in msg
                ]
                assert len(plotly_warnings) >= 1
                assert any("matplotlib" in msg.lower() for msg in plotly_warnings)

    def test_bokeh_fallback_to_matplotlib(self):
        """Falls back to matplotlib with warning when bokeh import fails."""
        from pylocuszoom.backends import _BACKENDS
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        # Clear bokeh from registry if present
        _BACKENDS.pop("bokeh", None)

        # Mock the bokeh imports to raise ImportError
        with patch.dict(
            "sys.modules",
            {
                "bokeh": None,
                "bokeh.plotting": None,
                "bokeh.models": None,
                "bokeh.layouts": None,
                "bokeh.io": None,
            },
        ):
            import sys

            sys.modules.pop("pylocuszoom.backends.bokeh_backend", None)

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                from importlib import reload

                import pylocuszoom.backends

                reload(pylocuszoom.backends)
                from pylocuszoom.backends import get_backend

                backend = get_backend("bokeh")

                # Should get matplotlib backend
                assert isinstance(backend, MatplotlibBackend)

                # Should have warned
                assert len(w) >= 1
                warning_messages = [str(warning.message) for warning in w]
                bokeh_warnings = [
                    msg
                    for msg in warning_messages
                    if "bokeh" in msg.lower() or "Bokeh" in msg
                ]
                assert len(bokeh_warnings) >= 1
                assert any("matplotlib" in msg.lower() for msg in bokeh_warnings)

    def test_fallback_warning_format(self):
        """Fallback warning mentions backend name and matplotlib."""
        from pylocuszoom.backends import _BACKENDS

        _BACKENDS.pop("plotly", None)

        with patch.dict(
            "sys.modules",
            {"plotly": None, "plotly.graph_objects": None, "plotly.subplots": None},
        ):
            import sys

            sys.modules.pop("pylocuszoom.backends.plotly_backend", None)

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")

                from importlib import reload

                import pylocuszoom.backends

                reload(pylocuszoom.backends)
                from pylocuszoom.backends import get_backend

                get_backend("plotly")

                # Find the specific warning about fallback
                fallback_warning = None
                for warning in w:
                    msg = str(warning.message).lower()
                    if "plotly" in msg and "matplotlib" in msg:
                        fallback_warning = str(warning.message)
                        break

                assert fallback_warning is not None, (
                    f"Expected fallback warning, got: {[str(w.message) for w in w]}"
                )
                # Check it mentions both the unavailable backend and fallback
                assert "plotly" in fallback_warning.lower()
                assert "matplotlib" in fallback_warning.lower()


class TestBackendCapabilities:
    """Tests that registered backends have expected capability properties."""

    def test_matplotlib_has_capabilities(self):
        """MatplotlibBackend has all capability properties."""
        from pylocuszoom.backends import get_backend

        backend = get_backend("matplotlib")

        assert hasattr(backend, "supports_snp_labels")
        assert hasattr(backend, "supports_hover")
        assert hasattr(backend, "supports_secondary_axis")

        # Matplotlib specific values
        assert backend.supports_snp_labels is True
        assert backend.supports_hover is False
        assert backend.supports_secondary_axis is True

    def test_plotly_has_capabilities(self):
        """PlotlyBackend has all capability properties."""
        pytest.importorskip("plotly")
        from pylocuszoom.backends import get_backend

        backend = get_backend("plotly")

        assert hasattr(backend, "supports_snp_labels")
        assert hasattr(backend, "supports_hover")
        assert hasattr(backend, "supports_secondary_axis")

        # Plotly specific values
        assert backend.supports_snp_labels is False
        assert backend.supports_hover is True
        assert backend.supports_secondary_axis is True

    def test_bokeh_has_capabilities(self):
        """BokehBackend has all capability properties."""
        pytest.importorskip("bokeh")
        from pylocuszoom.backends import get_backend

        backend = get_backend("bokeh")

        assert hasattr(backend, "supports_snp_labels")
        assert hasattr(backend, "supports_hover")
        assert hasattr(backend, "supports_secondary_axis")

        # Bokeh specific values
        assert backend.supports_snp_labels is False
        assert backend.supports_hover is True
        assert backend.supports_secondary_axis is True

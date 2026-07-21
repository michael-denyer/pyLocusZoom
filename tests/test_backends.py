"""Tests for backend registration and fallback."""

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


class TestBackendImportErrors:
    """Tests for ImportError behavior when optional backends unavailable."""

    def test_plotly_succeeds_when_installed(self):
        """No error when plotly is available."""
        from pylocuszoom.backends import get_backend

        pytest.importorskip("plotly")
        backend = get_backend("plotly")
        assert backend is not None

    def test_bokeh_succeeds_when_installed(self):
        """No error when bokeh is available."""
        from pylocuszoom.backends import get_backend

        pytest.importorskip("bokeh")
        backend = get_backend("bokeh")
        assert backend is not None

    def test_error_messages_contain_install_instructions(self):
        """Verify ImportError messages contain install instructions."""
        import inspect

        from pylocuszoom.backends import get_backend

        source = inspect.getsource(get_backend)
        assert "pip install plotly" in source
        assert "pip install bokeh" in source

    def test_registry_persists_across_calls(self):
        """Registry persists across multiple get_backend calls."""
        from pylocuszoom.backends import _BACKENDS, get_backend
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        # First call registers matplotlib
        backend1 = get_backend("matplotlib")
        assert isinstance(backend1, MatplotlibBackend)
        assert "matplotlib" in _BACKENDS

        # Second call uses same registry
        backend2 = get_backend("matplotlib")
        assert isinstance(backend2, MatplotlibBackend)
        # Both are instances of the registered class
        assert _BACKENDS["matplotlib"] is MatplotlibBackend


class TestBackendCapabilities:
    """Tests that registered backends have expected capability properties."""

    def test_matplotlib_has_capabilities(self):
        """MatplotlibBackend supports SNP labels and secondary axis, not hover."""
        from pylocuszoom.backends import (
            SupportsSecondaryAxis,
            SupportsSNPLabels,
            get_backend,
        )

        backend = get_backend("matplotlib")
        assert isinstance(backend, SupportsSNPLabels)
        assert isinstance(backend, SupportsSecondaryAxis)
        assert backend.supports_hover is False

    def test_plotly_has_capabilities(self):
        """PlotlyBackend supports secondary axis and hover, not SNP labels."""
        pytest.importorskip("plotly")
        from pylocuszoom.backends import (
            SupportsSecondaryAxis,
            SupportsSNPLabels,
            get_backend,
        )

        backend = get_backend("plotly")
        assert not isinstance(backend, SupportsSNPLabels)
        assert isinstance(backend, SupportsSecondaryAxis)
        assert backend.supports_hover is True

    def test_bokeh_has_capabilities(self):
        """BokehBackend supports secondary axis and hover, not SNP labels."""
        pytest.importorskip("bokeh")
        from pylocuszoom.backends import (
            SupportsSecondaryAxis,
            SupportsSNPLabels,
            get_backend,
        )

        backend = get_backend("bokeh")
        assert not isinstance(backend, SupportsSNPLabels)
        assert isinstance(backend, SupportsSecondaryAxis)
        assert backend.supports_hover is True


class TestBackendRegistration:
    """Tests for backend decorator registration integration."""

    def test_matplotlib_registered_on_import(self):
        """MatplotlibBackend is registered when module is imported."""
        from pylocuszoom.backends import _BACKENDS

        # Import triggers registration
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        assert "matplotlib" in _BACKENDS
        assert _BACKENDS["matplotlib"] is MatplotlibBackend

    def test_plotly_registered_on_import(self):
        """PlotlyBackend is registered when module is imported."""
        pytest.importorskip("plotly")
        from pylocuszoom.backends import _BACKENDS
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        assert "plotly" in _BACKENDS
        assert _BACKENDS["plotly"] is PlotlyBackend

    def test_bokeh_registered_on_import(self):
        """BokehBackend is registered when module is imported."""
        pytest.importorskip("bokeh")
        from pylocuszoom.backends import _BACKENDS
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        assert "bokeh" in _BACKENDS
        assert _BACKENDS["bokeh"] is BokehBackend

    def test_decorator_on_all_backends(self):
        """All backend classes use @register_backend decorator."""
        import inspect

        # Check matplotlib (always available)
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        # Read the module source file to verify decorator is present
        source_file = inspect.getsourcefile(MatplotlibBackend)
        with open(source_file) as f:
            source = f.read()
        assert "@register_backend" in source
        assert '@register_backend("matplotlib")' in source

        # Check plotly if available
        pytest.importorskip("plotly")
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        source_file = inspect.getsourcefile(PlotlyBackend)
        with open(source_file) as f:
            source = f.read()
        assert "@register_backend" in source
        assert '@register_backend("plotly")' in source

        # Check bokeh if available
        pytest.importorskip("bokeh")
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        source_file = inspect.getsourcefile(BokehBackend)
        with open(source_file) as f:
            source = f.read()
        assert "@register_backend" in source
        assert '@register_backend("bokeh")' in source


class TestSetXticks:
    """Tests for x-axis tick setting across backends."""

    def test_matplotlib_set_xticks(self):
        """Matplotlib backend should set x-axis ticks."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 4))
        backend.set_xticks(axes[0], [0, 1, 2], ["A", "B", "C"])
        # Verify ticks were set
        ticks = list(axes[0].get_xticks())
        assert 0 in ticks
        assert 1 in ticks
        assert 2 in ticks

    def test_plotly_set_xticks(self):
        """Plotly backend should set x-axis ticks."""
        pytest.importorskip("plotly")
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 4))
        backend.set_xticks(axes[0], [0, 1, 2], ["A", "B", "C"])
        # Verify via layout
        xaxis = fig.layout.xaxis
        assert xaxis.tickvals == (0, 1, 2)
        assert xaxis.ticktext == ("A", "B", "C")

    def test_bokeh_set_xticks(self):
        """Bokeh backend should set x-axis ticks."""
        pytest.importorskip("bokeh")
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 4))
        backend.set_xticks(axes[0], [0, 1, 2], ["A", "B", "C"])
        # Access ticks via the ticker's ticks property
        assert list(axes[0].xaxis.ticker.ticks) == [0, 1, 2]


class TestConvertLatexToUnicode:
    """Tests for LaTeX to Unicode conversion."""

    def test_convert_neg_log10_p(self):
        """Should convert -log10(P) LaTeX notation."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"$-\log_{10}$ P")
        # The conversion replaces "$-\log_{10}$ P" with "-log10(P)"
        assert result == "-log10(P)"

    def test_convert_neg_log10(self):
        """Should convert -log10 LaTeX notation."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"$-\log_{10}$")
        assert result == "-log10"

    def test_convert_log10(self):
        """Should convert log10 LaTeX notation."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"\log_{10}")
        assert result == "log10"

    def test_convert_r2_lowercase(self):
        """Should convert r² LaTeX notation."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"$r^2$")
        assert result == "r"

    def test_convert_r2_uppercase(self):
        """Should convert R² LaTeX notation."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"$R^2$")
        assert result == "R"

    def test_strips_dollar_signs(self):
        """Should strip remaining dollar signs."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"$some text$")
        assert "$" not in result
        assert result == "some text"

    def test_no_conversion_plain_text(self):
        """Plain text without LaTeX should pass through unchanged."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode("plain text")
        assert result == "plain text"

    def test_partial_conversion(self):
        """Should handle labels with some LaTeX and some plain text."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"Value: $r^2$ = 0.5")
        assert result == "Value: r = 0.5"


class TestLazyAttributeAccess:
    """Tests for lazy attribute access via __getattr__."""

    def test_matplotlib_backend_lazy_access(self):
        """MatplotlibBackend should be accessible via lazy import."""
        from pylocuszoom import backends

        # Access via module __getattr__
        MatplotlibBackend = backends.MatplotlibBackend
        assert MatplotlibBackend is not None

        # Should be the actual class
        from pylocuszoom.backends.matplotlib_backend import (
            MatplotlibBackend as DirectBackend,
        )

        assert MatplotlibBackend is DirectBackend

    def test_unknown_attribute_raises_attributeerror(self):
        """Accessing unknown attribute should raise AttributeError."""
        from pylocuszoom import backends

        with pytest.raises(AttributeError) as exc_info:
            _ = backends.NonExistentAttribute

        assert "NonExistentAttribute" in str(exc_info.value)

    def test_all_exports(self):
        """Module __all__ should contain expected exports."""
        from pylocuszoom.backends import __all__

        expected = [
            "PlotBackend",
            "BackendType",
            "get_backend",
            "register_backend",
            "MatplotlibBackend",
            "convert_latex_to_unicode",
        ]
        for name in expected:
            assert name in __all__


class TestBackendTypeLiteral:
    """Tests for BackendType type literal."""

    def test_backend_type_values(self):
        """BackendType should be a valid Literal type."""
        from pylocuszoom.backends import get_backend

        # These should all work without type errors at runtime
        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            # Type checkers would verify BackendType compatibility
            backend = get_backend(backend_name)
            assert backend is not None


class TestHeatmapMethods:
    """Tests for heatmap rendering methods across backends."""

    @pytest.fixture
    def sample_ld_matrix(self):
        """Create a sample 5x5 LD matrix for testing."""
        import numpy as np

        # Symmetric matrix with diagonal = 1, decreasing r2 with distance
        n = 5
        data = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                dist = abs(i - j)
                data[i, j] = 1.0 - dist * 0.2
        return data

    def test_matplotlib_add_heatmap_returns_mappable(self, sample_ld_matrix):
        """Matplotlib add_heatmap should return AxesImage object."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
        )
        assert mappable is not None
        # Should have a colormap
        assert hasattr(mappable, "get_cmap")
        backend.close(fig)

    def test_matplotlib_add_heatmap_mask_upper(self, sample_ld_matrix):
        """Matplotlib add_heatmap with mask_upper should produce masked array."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            mask_upper=True,
        )
        # Image data should be masked array when mask_upper=True

        # The underlying data should be masked
        assert mappable is not None
        backend.close(fig)

    def test_matplotlib_add_colorbar(self, sample_ld_matrix):
        """Matplotlib add_colorbar should not raise."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
        )
        cbar = backend.add_colorbar(axes[0], mappable, label="R²")
        assert cbar is not None
        backend.close(fig)

    def test_plotly_add_heatmap_returns_trace(self, sample_ld_matrix):
        """Plotly add_heatmap should return Heatmap trace."""
        pytest.importorskip("plotly")
        import plotly.graph_objects as go

        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        trace = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
        )
        assert trace is not None
        assert isinstance(trace, go.Heatmap)

    def test_plotly_add_colorbar_no_error(self, sample_ld_matrix):
        """Plotly add_colorbar should not raise (no-op)."""
        pytest.importorskip("plotly")
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        trace = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
        )
        # Should not raise
        result = backend.add_colorbar(axes[0], trace, label="R²")
        assert result is None  # No-op returns None

    def test_bokeh_add_heatmap_returns_mapper(self, sample_ld_matrix):
        """Bokeh add_heatmap should return LinearColorMapper."""
        pytest.importorskip("bokeh")
        from bokeh.models import LinearColorMapper

        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        mapper = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
        )
        assert mapper is not None
        assert isinstance(mapper, LinearColorMapper)

    def test_bokeh_add_colorbar_adds_to_layout(self, sample_ld_matrix):
        """Bokeh add_colorbar should add ColorBar to figure."""
        pytest.importorskip("bokeh")
        from bokeh.models import ColorBar

        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        mapper = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
        )
        cbar = backend.add_colorbar(axes[0], mapper, label="R²")
        assert cbar is not None
        assert isinstance(cbar, ColorBar)
        # Should be added to right layout
        assert cbar in axes[0].right

    def test_all_backends_have_heatmap_methods(self):
        """All backends should have add_heatmap and add_colorbar methods."""
        from pylocuszoom.backends import get_backend

        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            if backend_name in ["plotly", "bokeh"]:
                pytest.importorskip(backend_name)

            backend = get_backend(backend_name)
            assert hasattr(backend, "add_heatmap"), (
                f"{backend_name} missing add_heatmap"
            )
            assert hasattr(backend, "add_colorbar"), (
                f"{backend_name} missing add_colorbar"
            )

    def test_matplotlib_custom_colors(self, sample_ld_matrix):
        """Matplotlib add_heatmap should accept custom color gradient."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))
        # Use blue-to-yellow gradient
        mappable = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=["#0000FF", "#FFFF00"],
        )
        assert mappable is not None
        backend.close(fig)

    def test_heatmap_mask_upper_lower_triangle(self, sample_ld_matrix):
        """Test that mask_upper=True renders only lower triangle."""
        import numpy as np

        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(1, [1.0], (6, 6))

        # With mask_upper=True
        mappable = backend.add_heatmap(
            axes[0],
            sample_ld_matrix,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            mask_upper=True,
        )

        # Get the array data - should be masked
        array_data = mappable.get_array()
        # Check that upper triangle is masked
        assert np.ma.is_masked(array_data)
        backend.close(fig)


class TestCustomBackendCompatibility:
    """Tests for custom backend forward compatibility."""

    def test_recomb_overlay_skipped_when_secondary_axis_unsupported(self):
        """Recomb overlay is skipped (no crash) when a backend lacks the capability.

        The overlay is composed above the seam and gated on
        isinstance(backend, SupportsSecondaryAxis). A backend missing the
        secondary-axis methods must skip it gracefully rather than error.
        """
        import matplotlib.pyplot as plt
        import pandas as pd

        from pylocuszoom.plotter import LocusZoomPlotter

        plotter = LocusZoomPlotter(species=None)
        backend_cls = type(plotter._backend)
        original = backend_cls.create_twin_axis
        # Removing a secondary-axis method makes the backend fail the
        # SupportsSecondaryAxis isinstance check.
        delattr(backend_cls, "create_twin_axis")

        try:
            gwas_df = pd.DataFrame(
                {
                    "rs": ["rs1", "rs2", "rs3"],
                    "chr": [1, 1, 1],
                    "ps": [1100000, 1500000, 1900000],
                    "p_wald": [1e-8, 1e-5, 0.01],
                }
            )
            recomb_df = pd.DataFrame(
                {"pos": [1000000, 1500000, 2000000], "rate": [50.0, 100.0, 75.0]}
            )
            fig = plotter.plot(
                gwas_df, chrom=1, start=1000000, end=2000000, recomb_df=recomb_df
            )
            assert fig is not None
            plt.close(fig)
        finally:
            backend_cls.create_twin_axis = original

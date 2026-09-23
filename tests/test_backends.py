"""Tests for backend registration and fallback."""

import sys

import numpy as np
import pytest

from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.colors import LD_HEATMAP_COLORS
from tests.figure_probes import PROBES


class TestRegisterBackend:
    """Tests for the @register_backend decorator."""

    def test_registered_backend_resolves_through_get_backend(self, monkeypatch):
        """A name registered by the decorator is resolvable by that name."""
        from pylocuszoom.backends import _BACKENDS, get_backend, register_backend

        monkeypatch.delitem(_BACKENDS, "test_dummy", raising=False)

        @register_backend("test_dummy")
        class DummyBackend:
            pass

        assert isinstance(get_backend("test_dummy"), DummyBackend)

    def test_register_backend_returns_class_unchanged(self, monkeypatch):
        """Decorator returns the class unchanged."""
        from pylocuszoom.backends import _BACKENDS, register_backend

        monkeypatch.delitem(_BACKENDS, "test_unchanged", raising=False)

        @register_backend("test_unchanged")
        class OriginalBackend:
            def method(self):
                return "original"

        assert OriginalBackend().method() == "original"


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


class TestBackendImportErrors:
    """Tests for ImportError behavior when optional backends unavailable."""

    @pytest.mark.parametrize("name", ["plotly", "bokeh"])
    def test_missing_optional_backend_names_its_install_command(
        self, monkeypatch, name
    ):
        """An uninstalled optional backend raises ImportError naming the pip command."""
        from pylocuszoom import backends

        for cached in [m for m in sys.modules if m.split(".")[0] == name]:
            monkeypatch.delitem(sys.modules, cached)
        monkeypatch.setitem(sys.modules, name, None)
        monkeypatch.delitem(
            sys.modules, f"pylocuszoom.backends.{name}_backend", raising=False
        )
        monkeypatch.delitem(backends._BACKENDS, name, raising=False)

        with pytest.raises(ImportError, match=f"pip install {name}"):
            backends.get_backend(name)


class TestBackendCapabilities:
    """Tests that registered backends have expected capability properties."""

    @pytest.mark.parametrize(
        ("backend_name", "labels", "hover"),
        [("matplotlib", True, False), ("plotly", False, True), ("bokeh", False, True)],
    )
    def test_backend_declares_its_capabilities(self, backend_name, labels, hover):
        """Only matplotlib labels SNPs; only the interactive backends hover."""
        from pylocuszoom.backends import SupportsSNPLabels, get_backend

        backend = get_backend(backend_name)
        assert isinstance(backend, SupportsSNPLabels) is labels
        assert backend.supports_hover is hover


class TestBackendRegistration:
    """Every built-in backend reaches the registry through its decorator."""

    @staticmethod
    def _builtin_classes():
        from pylocuszoom.backends.bokeh_backend import BokehBackend
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        return {
            "matplotlib": MatplotlibBackend,
            "plotly": PlotlyBackend,
            "bokeh": BokehBackend,
        }

    @pytest.mark.parametrize("name", BUILTIN_BACKENDS)
    def test_builtin_name_resolves_to_its_backend_class(self, name):
        """Importing a backend module makes its name resolvable to that class."""
        from pylocuszoom.backends import get_backend

        assert isinstance(get_backend(name), self._builtin_classes()[name])

    def test_builtin_backends_covers_every_registered_name(self):
        """BUILTIN_BACKENDS names exactly the backends that ship with the package."""
        assert set(BUILTIN_BACKENDS) == set(self._builtin_classes())


class TestSetXticks:
    """Tests for x-axis tick setting across backends."""

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_set_xticks_places_positions_and_labels(self, backend_name):
        """set_xticks puts each label at its position."""
        from pylocuszoom.backends import get_backend

        backend = get_backend(backend_name)
        fig, axes = backend.create_figure([1.0], (6, 4))
        backend.set_xticks(axes[0], [0, 1, 2], ["A", "B", "C"])

        assert PROBES[backend_name].xticks(fig) == ([0, 1, 2], ["A", "B", "C"])


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
        assert result == "r²"

    def test_convert_r2_uppercase(self):
        """Should convert R² LaTeX notation."""
        from pylocuszoom.backends import convert_latex_to_unicode

        result = convert_latex_to_unicode(r"$R^2$")
        assert result == "R²"

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
        assert result == "Value: r² = 0.5"


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


class TestHeatmapMethods:
    """Tests for heatmap rendering methods across backends."""

    @pytest.fixture
    def ld_matrix_array(self):
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

    def test_matplotlib_add_heatmap_returns_mappable(self, ld_matrix_array):
        """Matplotlib add_heatmap returns a mesh usable by a colorbar."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        assert mappable is not None
        # Should have a colormap
        assert hasattr(mappable, "get_cmap")

    def test_matplotlib_add_heatmap_lower_triangle(self, ld_matrix_array):
        """Draw a lower-triangle matrix with its upper triangle masked out."""
        from pylocuszoom.backends.composition import lower_triangle
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            lower_triangle(ld_matrix_array),
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )

        drawn = np.ma.getmaskarray(mappable.get_array()).reshape(5, 5)
        assert drawn.tolist() == np.triu(np.ones((5, 5), dtype=bool), k=1).tolist()

    def test_matplotlib_add_colorbar(self, ld_matrix_array):
        """Matplotlib add_colorbar attaches a labelled scale to the figure."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        backend.add_colorbar(axes[0], mappable, label="R²")

        assert [a.get_ylabel() for a in fig.axes if a is not axes[0]] == ["R²"]

    def test_plotly_add_heatmap_returns_trace(self, ld_matrix_array):
        """Plotly add_heatmap should return Heatmap trace."""
        import plotly.graph_objects as go

        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        trace = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        assert trace is not None
        assert isinstance(trace, go.Heatmap)

    def test_plotly_add_colorbar_enables_the_trace_scale(self, ld_matrix_array):
        """Plotly's colorbar is the trace's own scale, off until asked for."""
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        trace = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        assert trace.showscale is False

        backend.add_colorbar(axes[0], trace, label="D'")

        assert trace.showscale is True
        assert trace.colorbar.title.text == "D'"

    def test_plotly_add_colorbar_honours_orientation(self, ld_matrix_array):
        """Horizontal orientation maps to Plotly's 'h'."""
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        trace = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        backend.add_colorbar(axes[0], trace, label="R²", orientation="horizontal")

        assert trace.colorbar.orientation == "h"

    def test_bokeh_add_heatmap_returns_mapper(self, ld_matrix_array):
        """Bokeh add_heatmap should return LinearColorMapper."""
        from bokeh.models import LinearColorMapper

        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        mapper = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        assert mapper is not None
        assert isinstance(mapper, LinearColorMapper)

    def test_bokeh_add_colorbar_adds_to_layout(self, ld_matrix_array):
        """Bokeh add_colorbar should add ColorBar to figure."""
        from bokeh.models import ColorBar

        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        mapper = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )
        backend.add_colorbar(axes[0], mapper, label="R²")

        assert [type(m) for m in axes[0].right] == [ColorBar]

    def test_matplotlib_custom_colors(self, ld_matrix_array):
        """Build the heatmap colormap from the caller's own gradient stops."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))
        mappable = backend.add_heatmap(
            axes[0],
            ld_matrix_array,
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=["#0000FF", "#FFFF00"],
        )

        cmap = mappable.get_cmap()
        assert cmap(0.0) == (0.0, 0.0, 1.0, 1.0)
        assert cmap(1.0) == (1.0, 1.0, 0.0, 1.0)

    def test_heatmap_lower_triangle_masks_upper(self, ld_matrix_array):
        """A lower_triangle matrix should reach matplotlib still masked."""
        import numpy as np

        from pylocuszoom.backends.composition import lower_triangle
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure([1.0], (6, 6))

        mappable = backend.add_heatmap(
            axes[0],
            lower_triangle(ld_matrix_array),
            x_coords=list(range(5)),
            y_coords=list(range(5)),
            cmap_colors=LD_HEATMAP_COLORS,
        )

        # Get the array data - should be masked
        array_data = mappable.get_array()
        # Check that upper triangle is masked
        assert np.ma.is_masked(array_data)


class TestCustomBackendCompatibility:
    """Tests for custom backend forward compatibility."""

    def test_backend_missing_snp_labels_lacks_the_one_optional_capability(self):
        """SNP labels are the only capability a backend can decline."""
        from pylocuszoom.backends import SupportsSNPLabels

        class MinimalBackend:
            pass

        class LabellingBackend:
            def add_snp_labels(self, *args, **kwargs): ...

        assert not isinstance(MinimalBackend(), SupportsSNPLabels)
        assert isinstance(LabellingBackend(), SupportsSNPLabels)


class TestLegendPlacement:
    """Every backend honours the neutral add_legend contract."""

    @staticmethod
    def _entries():
        from pylocuszoom.backends.composition import LegendEntry

        return [
            LegendEntry("Lead SNP", "#FF0000", marker="D", edgecolor="#00FF00"),
            LegendEntry("0.8 - 1.0", "#0000FF"),
        ]

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_honours_loc(self, backend_name):
        """The legend is anchored in the corner the matplotlib loc names."""
        from pylocuszoom.backends import get_backend

        backend = get_backend(backend_name)
        fig, axes = backend.create_figure([1.0], (6, 4))
        backend.add_legend(axes[0], self._entries(), loc="lower left")

        assert PROBES[backend_name].legend_corner(fig) == "lower left"

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_honours_edgecolor(self, backend_name):
        """A swatch edge takes its entry's edgecolor, black when it has none."""
        from pylocuszoom.backends import get_backend

        backend = get_backend(backend_name)
        fig, axes = backend.create_figure([1.0], (6, 4))
        backend.add_legend(axes[0], self._entries(), loc="upper right")

        assert PROBES[backend_name].legend_edgecolors(fig) == {
            "Lead SNP": "#00ff00",
            "0.8 - 1.0": "#000000",
        }

    def test_unknown_loc_falls_back_without_raising(self):
        """An unmapped loc degrades to the default corner, it does not raise."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        plotly_backend = PlotlyBackend()
        plotly_fig, plotly_axes = plotly_backend.create_figure([1.0], (6, 4))
        plotly_backend.add_legend(plotly_axes[0], self._entries(), loc="nonsense")
        assert plotly_fig.layout.legend.xanchor == "right"

        bokeh_backend = BokehBackend()
        _, bokeh_axes = bokeh_backend.create_figure([1.0], (6, 4))
        bokeh_backend.add_legend(bokeh_axes[0], self._entries(), loc="nonsense")


class TestLegendTitleMathtext:
    """The LD legend title is mathtext, rendered natively per backend."""

    def test_matplotlib_keeps_mathtext(self):
        """Matplotlib receives the raw mathtext so it renders an italic r²."""
        from pylocuszoom.backends.composition import (
            LD_LEGEND_TITLE,
            ld_legend_entries,
        )
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        assert LD_LEGEND_TITLE == r"$r^2$"

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure([1.0], (6, 4))
        backend.add_legend(
            axes[0], ld_legend_entries(), loc="upper right", title=LD_LEGEND_TITLE
        )
        assert axes[0].get_legend().get_title().get_text() == r"$r^2$"

    def test_interactive_backends_show_unicode(self):
        """Plotly and Bokeh convert the mathtext to a plain unicode r²."""
        from bokeh.models import Legend

        from pylocuszoom.backends.bokeh_backend import BokehBackend
        from pylocuszoom.backends.composition import (
            LD_LEGEND_TITLE,
            ld_legend_entries,
        )
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        plotly_backend = PlotlyBackend()
        plotly_fig, plotly_axes = plotly_backend.create_figure([1.0], (6, 4))
        plotly_backend.add_legend(
            plotly_axes[0], ld_legend_entries(), title=LD_LEGEND_TITLE
        )
        assert plotly_fig.layout.legend.title.text == "r²"

        bokeh_backend = BokehBackend()
        _, bokeh_axes = bokeh_backend.create_figure([1.0], (6, 4))
        bokeh_backend.add_legend(
            bokeh_axes[0], ld_legend_entries(), title=LD_LEGEND_TITLE
        )
        assert list(bokeh_axes[0].select(Legend))[0].title == "r²"


@pytest.mark.parametrize("orientation", ["vertical", "horizontal"])
def test_matplotlib_standalone_colorbar_survives_final_layout(orientation):
    from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

    backend = MatplotlibBackend()
    fig, axes = backend.create_figure([1.0], (6, 6))
    mappable = backend.add_heatmap(axes[0], np.eye(2), [0, 1], [0, 1], ["white", "red"])
    backend.add_colorbar(axes[0], mappable, label="Scale", orientation=orientation)
    backend.finalize_layout(fig)
    fig.canvas.draw()
    colorbar = fig.axes[1]
    bounds = colorbar.get_position()
    if orientation == "vertical":
        assert colorbar.get_ylabel() == "Scale"
        assert bounds.height > bounds.width
        assert bounds.x0 > axes[0].get_position().x1
    else:
        assert colorbar.get_xlabel() == "Scale"
        assert bounds.width > bounds.height
        assert bounds.y1 < axes[0].get_position().y0

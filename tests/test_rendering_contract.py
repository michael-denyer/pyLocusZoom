"""Contract tests for the semantic Manhattan/QQ renderer."""

from types import SimpleNamespace

import pandas as pd
import pytest

from pylocuszoom._miami_renderer import MiamiRenderer
from pylocuszoom._rendering import ManhattanQQRenderer
from pylocuszoom.backends import get_backend
from pylocuszoom.manhattan import prepare_manhattan_data
from pylocuszoom.qq import prepare_qq_data


class RecordingBackend:
    """Small primitive adapter used to test the renderer's semantic seam."""

    supports_hover = False

    def __init__(self):
        self.calls = []

    def _record(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))
        return SimpleNamespace()

    def create_figure(self, **kwargs):
        self._record("create_figure", **kwargs)
        return SimpleNamespace(), [SimpleNamespace() for _ in range(kwargs["n_panels"])]

    def create_figure_grid(self, **kwargs):
        self._record("create_figure_grid", **kwargs)
        n_axes = kwargs["n_rows"] * kwargs["n_cols"]
        return SimpleNamespace(), [SimpleNamespace() for _ in range(n_axes)]

    def finalize_layout(self, *args, **kwargs):
        self._record("finalize_layout", *args, **kwargs)

    def scatter(self, *args, **kwargs):
        return self._record("scatter", *args, **kwargs)

    def fill_between(self, *args, **kwargs):
        return self._record("fill_between", *args, **kwargs)

    def line(self, *args, **kwargs):
        return self._record("line", *args, **kwargs)

    def axhline(self, *args, **kwargs):
        return self._record("axhline", *args, **kwargs)

    def set_xlim(self, *args, **kwargs):
        self._record("set_xlim", *args, **kwargs)

    def set_ylim(self, *args, **kwargs):
        self._record("set_ylim", *args, **kwargs)

    def set_xticks(self, *args, **kwargs):
        self._record("set_xticks", *args, **kwargs)

    def set_xlabel(self, *args, **kwargs):
        self._record("set_xlabel", *args, **kwargs)

    def set_ylabel(self, *args, **kwargs):
        self._record("set_ylabel", *args, **kwargs)

    def set_title(self, *args, **kwargs):
        self._record("set_title", *args, **kwargs)

    def set_suptitle(self, *args, **kwargs):
        self._record("set_suptitle", *args, **kwargs)

    def hide_spines(self, *args, **kwargs):
        self._record("hide_spines", *args, **kwargs)

    def add_panel_label(self, *args, **kwargs):
        self._record("add_panel_label", *args, **kwargs)


@pytest.fixture
def prepared_data():
    df = pd.DataFrame(
        {
            "chrom": [1, 1, 2, 2],
            "pos": [1_000_000, 2_000_000, 1_000_000, 2_000_000],
            "p": [1e-8, 0.01, 1e-6, 0.5],
        }
    )
    return (
        prepare_manhattan_data(df, species="human"),
        prepare_qq_data(df),
    )


def test_renderer_owns_manhattan_panel_policy(prepared_data):
    manhattan_df, _ = prepared_data
    backend = RecordingBackend()

    ManhattanQQRenderer(backend).render_manhattan(
        manhattan_df,
        figsize=(12, 5),
        significance_threshold=5e-8,
        title="Contract Manhattan",
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert "scatter" in names
    assert "axhline" in names
    assert names[-1] == "finalize_layout"
    assert "set_title" in names
    assert "set_xticks" in names


def test_renderer_owns_qq_panel_policy(prepared_data):
    _, qq_df = prepared_data
    backend = RecordingBackend()

    ManhattanQQRenderer(backend).render_qq(
        qq_df,
        figsize=(6, 6),
        show_confidence_band=True,
        show_lambda=True,
        title=None,
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert "fill_between" in names
    assert "line" in names
    assert "scatter" in names
    assert names[-1] == "finalize_layout"


@pytest.mark.parametrize("backend_name", ["matplotlib", "plotly", "bokeh"])
def test_same_prepared_intent_renders_on_each_builtin_backend(
    backend_name, prepared_data
):
    """All built-in adapters accept the same prepared rendering intent."""
    if backend_name != "matplotlib":
        pytest.importorskip(backend_name)

    manhattan_df, qq_df = prepared_data
    backend = get_backend(backend_name)
    renderer = ManhattanQQRenderer(backend)
    figures = []
    try:
        figures.append(
            renderer.render_manhattan(
                manhattan_df,
                figsize=(12, 5),
                significance_threshold=5e-8,
                title="Contract Manhattan",
            )
        )
        figures.append(
            renderer.render_qq(
                qq_df,
                figsize=(6, 6),
                show_confidence_band=True,
                show_lambda=True,
                title=None,
            )
        )
        assert all(figure is not None for figure in figures)
    finally:
        for figure in figures:
            backend.close(figure)


def test_miami_region_highlight_is_optional_for_legacy_backends(prepared_data):
    """A pre-existing custom backend need not implement the new capability."""
    manhattan_df, _ = prepared_data
    backend = RecordingBackend()

    figure = MiamiRenderer(backend).render(
        manhattan_df,
        manhattan_df,
        pos_col="pos",
        p_col="p",
        rs_col=None,
        top_threshold=None,
        bottom_threshold=None,
        top_label=None,
        bottom_label=None,
        top_snp_annotations=None,
        bottom_snp_annotations=None,
        highlight_regions=[("1", 1_000_000, 2_000_000)],
        highlight_color="yellow",
        highlight_alpha=0.3,
        figsize=(12, 8),
        title=None,
    )

    assert figure is not None

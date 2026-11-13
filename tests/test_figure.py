"""render_figure issues exactly the figure-level primitives a plan names."""

import pandas as pd
import pytest

from pylocuszoom._figure import FigurePlan, RegionHighlight, render_figure
from pylocuszoom.backends import BUILTIN_BACKENDS, get_backend
from tests.test_rendering_contract import RecordingBackend


class _Panel:
    """A panel that records the axis it was drawn on."""

    def __init__(self):
        self.drawn_on = []

    def draw(self, backend, ax):
        self.drawn_on.append(ax)
        backend.scatter(ax, [], [])


def _calls(backend, name):
    return [(args, kwargs) for called, args, kwargs in backend.calls if called == name]


def test_column_plan_creates_one_axis_per_panel_and_draws_in_order():
    backend = RecordingBackend()
    panels = [_Panel(), _Panel(), _Panel()]

    render_figure(
        backend,
        FigurePlan(panels=panels, figsize=(8.0, 6.0), height_ratios=[3.0, 1.0, 1.0]),
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert names[-1] == "finalize_layout"
    assert "create_figure_grid" not in names
    (_, kwargs), *_ = _calls(backend, "create_figure")
    assert kwargs == {
        "height_ratios": [3.0, 1.0, 1.0],
        "figsize": (8.0, 6.0),
        "sharex": True,
    }
    assert [len(panel.drawn_on) for panel in panels] == [1, 1, 1]
    assert len({id(panel.drawn_on[0]) for panel in panels}) == 3


def test_equal_rows_when_no_height_ratios_are_given():
    backend = RecordingBackend()

    render_figure(backend, FigurePlan(panels=[_Panel(), _Panel()], figsize=(1, 1)))

    (_, kwargs), *_ = _calls(backend, "create_figure")
    assert kwargs["height_ratios"] == [1.0, 1.0]


def test_grid_plan_uses_create_figure_grid():
    backend = RecordingBackend()
    panels = [_Panel() for _ in range(4)]

    render_figure(
        backend,
        FigurePlan(panels=panels, figsize=(14.0, 8.0), n_cols=2, width_ratios=[2.5, 1]),
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure_grid"
    assert "create_figure" not in names
    (_, kwargs), *_ = _calls(backend, "create_figure_grid")
    assert kwargs == {
        "n_rows": 2,
        "n_cols": 2,
        "width_ratios": [2.5, 1],
        "height_ratios": None,
        "figsize": (14.0, 8.0),
    }
    assert [len(panel.drawn_on) for panel in panels] == [1, 1, 1, 1]


def test_figure_level_policy_reaches_the_backend():
    backend = RecordingBackend()

    render_figure(
        backend,
        FigurePlan(
            panels=[_Panel(), _Panel()],
            figsize=(12.0, 8.0),
            xlabel="Chromosome 1 (Mb)",
            mb_xaxis=True,
            highlights=[RegionHighlight(10.0, 20.0, "yellow", 0.3)],
            suptitle="Title",
            first_panel_title="Top panel",
            top=0.9,
            hspace=0.05,
        ),
    )

    _, axes = backend.create_figure(height_ratios=[1.0, 1.0], figsize=(1, 1))
    ((xlabel_args, _),) = _calls(backend, "set_xlabel")
    assert xlabel_args[1] == "Chromosome 1 (Mb)"
    names = [name for name, _, _ in backend.calls]
    assert names.count("format_xaxis_mb") == 2
    assert names.index("set_xlabel") < names.index("format_xaxis_mb")
    ((highlight_args, highlight_kwargs),) = _calls(backend, "add_region_highlight")
    assert len(highlight_args[0]) == 2
    assert highlight_args[1:] == (10.0, 20.0)
    assert highlight_kwargs == {"color": "yellow", "alpha": 0.3}
    ((suptitle_args, suptitle_kwargs),) = _calls(backend, "set_suptitle")
    assert suptitle_args[1] == "Title" and suptitle_kwargs == {
        "fontsize": 14,
        "fontweight": "bold",
    }
    ((title_args, title_kwargs),) = _calls(backend, "set_title")
    assert title_args[1] == "Top panel" and title_kwargs == {
        "fontsize": 14,
        "fontweight": "bold",
    }
    assert names.index("scatter", 1) < names.index("set_title")
    ((_, layout_kwargs),) = _calls(backend, "finalize_layout")
    assert layout_kwargs == {"top": 0.9, "hspace": 0.05}


def test_no_title_means_no_suptitle_call():
    backend = RecordingBackend()

    render_figure(backend, FigurePlan(panels=[_Panel()], figsize=(1, 1)))

    names = [name for name, _, _ in backend.calls]
    assert "set_suptitle" not in names
    assert "set_title" not in names
    assert "set_xlabel" not in names
    assert "format_xaxis_mb" not in names
    assert "add_region_highlight" not in names


def test_empty_plan_is_rejected():
    with pytest.raises(ValueError, match="at least one panel"):
        render_figure(RecordingBackend(), FigurePlan(panels=[], figsize=(1, 1)))


def test_figure_level_calls_follow_the_draw_in_a_fixed_order():
    """Plotly writes layout keys in call order, so the order is the export."""
    backend = RecordingBackend()

    render_figure(
        backend,
        FigurePlan(
            panels=[_Panel(), _Panel()],
            figsize=(12.0, 8.0),
            xlabel="Chromosome 1 (Mb)",
            mb_xaxis=True,
            highlights=[RegionHighlight(10.0, 20.0, "yellow", 0.3)],
            suptitle="Title",
            first_panel_title="Top panel",
            footer="n = 1,234",
        ),
    )

    assert [name for name, _, _ in backend.calls] == [
        "create_figure",
        "scatter",
        "scatter",
        "set_xlabel",
        "format_xaxis_mb",
        "format_xaxis_mb",
        "add_region_highlight",
        "set_title",
        "set_suptitle",
        "finalize_layout",
        "set_footer",
    ]
    ((footer_args, _),) = _calls(backend, "set_footer")
    assert footer_args[1] == "n = 1,234"


class _DotPanel:
    """A panel that draws one point, so every real backend has something to lay out."""

    def draw(self, backend, ax):
        backend.scatter(ax, pd.Series([1.0]), pd.Series([1.0]), colors="black")


def _bottom_margin(backend_name, fig):
    """The space a backend leaves under its panels, in its own units."""
    if backend_name == "matplotlib":
        return fig.subplotpars.bottom
    if backend_name == "plotly":
        return fig.layout.margin.b
    raise AssertionError(f"no bottom margin for {backend_name}")


@pytest.mark.parametrize("backend_name", ["matplotlib", "plotly"])
def test_footer_room_survives_the_layout(backend_name):
    """set_footer grows the margin finalize_layout set; before it, the room is lost."""
    backend = get_backend(backend_name)
    figsize = (8.0, 6.0)

    def render(footer):
        plan = FigurePlan(panels=[_DotPanel()], figsize=figsize, footer=footer)
        return _bottom_margin(backend_name, render_figure(backend, plan))

    grown = render("n = 1,234") - render(None)

    # One 10pt line: 2 * fontsize, in pixels for plotly, in figure height for matplotlib.
    one_line = 20 if backend_name == "plotly" else 20 / 72 / figsize[1]
    assert grown == pytest.approx(one_line)


@pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
def test_render_writes_exactly_one_footer(backend_name):
    """set_footer is not idempotent on any backend, so render_figure calls it once."""
    fig = render_figure(
        get_backend(backend_name),
        FigurePlan(panels=[_DotPanel()], figsize=(8.0, 6.0), footer="n = 1,234"),
    )

    if backend_name == "matplotlib":
        texts = [text.get_text() for text in fig.texts]
    elif backend_name == "plotly":
        texts = [annotation.text for annotation in fig.layout.annotations]
    else:
        from bokeh.models import Div

        texts = [div.text for div in fig.select({"type": Div})]
    assert [text for text in texts if "n = 1,234" in text] == [texts[-1]], (
        f"expected one footer, got {texts}"
    )

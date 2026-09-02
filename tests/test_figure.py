"""render_figure issues exactly the figure-level primitives a plan names."""

import pytest

from pylocuszoom._figure import FigurePlan, RegionHighlight, render_figure
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
    assert suptitle_args[1] == "Title" and suptitle_kwargs == {"fontsize": 14}
    ((_, layout_kwargs),) = _calls(backend, "finalize_layout")
    assert layout_kwargs == {"top": 0.9, "hspace": 0.05}


def test_no_title_means_no_suptitle_call():
    backend = RecordingBackend()

    render_figure(backend, FigurePlan(panels=[_Panel()], figsize=(1, 1)))

    names = [name for name, _, _ in backend.calls]
    assert "set_suptitle" not in names
    assert "set_xlabel" not in names
    assert "format_xaxis_mb" not in names
    assert "add_region_highlight" not in names


def test_empty_plan_is_rejected():
    with pytest.raises(ValueError, match="at least one panel"):
        render_figure(RecordingBackend(), FigurePlan(panels=[], figsize=(1, 1)))

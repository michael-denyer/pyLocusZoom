"""Regional requests select one set of rows and resolve panel options once."""

from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pylocuszoom import ColumnConfig, DisplayConfig, LDConfig, LocusZoomPlotter

DISPLAY = DisplayConfig(show_recombination=False, snp_labels=False)


def points_and_limit(fig, backend):
    if backend == "matplotlib":
        return fig.axes[0].collections[0].get_offsets().tolist(), fig.axes[
            0
        ].get_ylim()[1]
    if backend == "plotly":
        return list(
            map(list, zip(fig.data[0].x, fig.data[0].y))
        ), fig.layout.yaxis.range[1]
    ax = fig.children[0]
    data = ax.renderers[0].data_source.data
    return list(map(list, zip(data["x"], data["y"]))), ax.y_range.end


@pytest.mark.parametrize("backend", ["matplotlib", "plotly", "bokeh"])
def test_region_selects_points_and_scale_before_drawing(backend):
    frame = pd.DataFrame(
        {
            "chr": [2, 1, 1],
            "pos": [150, 150, 9_000_000],
            "p_value": [1e-8, 0.01, 1e-200],
            "rs": ["wrong", "wanted", "far"],
        }
    )
    original = frame.copy(deep=True)
    plotter = LocusZoomPlotter(species=None, backend=backend, log_level=None)
    fig = plotter.plot(frame, chrom=1, start=100, end=200, display=DISPLAY)
    points, maximum = points_and_limit(fig, backend)
    assert points == [[150, 2.0]]
    assert maximum == pytest.approx(2.3)
    pd.testing.assert_frame_equal(frame, original)
    plt.close("all")


def test_ld_lead_comes_from_selected_chromosome():
    frame = pd.DataFrame(
        {
            "chr": [2, 1],
            "pos": [150, 150],
            "p_value": [1e-8, 0.01],
            "rs": ["wrong", "wanted"],
        }
    )
    leads = []

    def calculate(**kwargs):
        leads.append(kwargs["lead_snp"])
        return pd.DataFrame({"SNP": ["wanted"], "R2": [1.0]})

    with patch("pylocuszoom._ld_plotting.calculate_ld", side_effect=calculate):
        LocusZoomPlotter(species=None, log_level=None).plot(
            frame,
            chrom=1,
            start=100,
            end=200,
            display=DISPLAY,
            ld=LDConfig(lead_pos=150, ld_reference_file="reference"),
        )
    assert leads == ["wanted"]


@pytest.mark.parametrize("override, expected", [(None, 150), ([175], 175)])
def test_stacked_resolves_shared_lead_and_per_panel_override(override, expected):
    frame = pd.DataFrame({"pos": [150, 175], "p_value": [0.01, 1e-8]})
    fig = LocusZoomPlotter(species=None, backend="plotly", log_level=None).plot_stacked(
        [frame],
        chrom=1,
        start=100,
        end=200,
        display=DISPLAY,
        ld=LDConfig(lead_pos=150),
        lead_positions=override,
    )
    assert list(fig.data[-1].x) == [expected]


@pytest.mark.parametrize("lead", [0, -1])
def test_stacked_rejects_invalid_lead_positions(lead):
    frame = pd.DataFrame({"pos": [150], "p_value": [0.01]})
    with pytest.raises(ValueError, match="greater than or equal to 1"):
        LocusZoomPlotter(species=None, log_level=None).plot_stacked(
            [frame], chrom=1, start=100, end=200, display=DISPLAY, lead_positions=[lead]
        )


def test_stacked_resolves_columns_per_frame():
    canonical = pd.DataFrame({"pos": [150], "p_value": [0.01]})
    legacy = pd.DataFrame({"ps": [160], "p_wald": [0.001]})
    with pytest.warns(DeprecationWarning):
        fig = LocusZoomPlotter(
            species=None, backend="plotly", log_level=None
        ).plot_stacked(
            [canonical, legacy], chrom=1, start=100, end=200, display=DISPLAY
        )
    assert list(fig.data[0].x) == [150]
    assert list(fig.data[2].x) == [160]


def test_duplicate_index_does_not_make_lead_ambiguous():
    frame = pd.DataFrame({"pos": [150, 175], "p_value": [0.01, 1e-8]}, index=[0, 0])
    fig = LocusZoomPlotter(species=None, backend="plotly", log_level=None).plot(
        frame, chrom=1, start=100, end=200, display=DISPLAY, columns=ColumnConfig()
    )
    assert list(fig.data[-1].x) == [175]


@pytest.mark.parametrize("lead_positions", [None, [150]])
def test_duplicate_positions_use_one_strongest_lead_row(lead_positions):
    frame = pd.DataFrame(
        {
            "chr": [1, 1],
            "pos": [150, 150],
            "p_value": [0.1, 1e-8],
            "rs": ["weak", "strong"],
        }
    )
    observed_leads = []

    def calculate(**kwargs):
        observed_leads.append(kwargs["lead_snp"])
        return pd.DataFrame({"SNP": ["weak", "strong"], "R2": [0.2, 1.0]})

    with patch("pylocuszoom._ld_plotting.calculate_ld", side_effect=calculate):
        fig = LocusZoomPlotter(
            species=None, backend="plotly", log_level=None
        ).plot_stacked(
            [frame],
            chrom=1,
            start=100,
            end=200,
            display=DISPLAY,
            ld_reference_files=["reference"],
            lead_positions=lead_positions,
        )
    assert observed_leads == ["strong"]
    lead_trace = next(
        trace
        for trace in fig.data
        if trace.marker.symbol == "diamond" and trace.x[0] is not None
    )
    assert list(lead_trace.x) == [150]
    assert list(lead_trace.y) == [8.0]
    assert "strong" in str(lead_trace.customdata)


def test_requested_lead_missing_from_region_warns(warning_records):
    frame = pd.DataFrame({"pos": [150], "p_value": [0.1], "rs": ["other"]})
    fig = LocusZoomPlotter(species=None, backend="plotly", log_level=None).plot(
        frame,
        chrom=1,
        start=100,
        end=200,
        display=DISPLAY,
        ld=LDConfig(lead_pos=175, ld_reference_file="unused"),
    )
    assert list(fig.data[0].x) == [150]
    assert any("175 not found in region" in message for message in warning_records)


def test_same_position_nonlead_variant_is_excluded_from_lead_labels():
    frame = pd.DataFrame(
        {"pos": [150, 150], "p_value": [0.1, 1e-8], "rs": ["weak", "strong"]}
    )
    fig = LocusZoomPlotter(species=None, log_level=None).plot(
        frame,
        chrom=1,
        start=100,
        end=200,
        display=DisplayConfig(show_recombination=False, label_top_n=2),
    )
    assert [text.get_text() for text in fig.axes[0].texts] == ["strong"]

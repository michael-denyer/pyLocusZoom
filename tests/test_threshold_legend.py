"""The key to the threshold lines on Manhattan-family panels."""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom import GenomeWideStyle, ManhattanPlotter, MiamiPlotter
from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.backends.composition import (
    LegendEntry,
    threshold_label,
    threshold_legend_entries,
)
from tests.figure_probes import PROBES

SIGNIFICANCE = ("P = 5e-08", "#ff0000", "--")
NO_KEY = GenomeWideStyle(show_threshold_legend=False)


def _gwas(seed):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "chr": np.repeat([1, 2, 3, 4], 10),
            "pos": np.tile(np.arange(1, 11) * 1_000_000, 4),
            "p_value": np.concatenate(
                [rng.uniform(1e-3, 1, 10), [1e-9], rng.uniform(1e-3, 1, 29)]
            ),
        }
    )


@pytest.fixture
def threshold_gwas_df():
    return _gwas(7)


@pytest.fixture
def second_threshold_gwas_df():
    return _gwas(11)


class TestThresholdLabel:
    @pytest.mark.parametrize(
        ("threshold", "label"),
        [
            (5e-8, "P = 5e-08"),
            (1e-7, "P = 1e-07"),
            (1e-5, "P = 1e-05"),
            (0.05, "P = 5e-02"),
            (2.5e-6, "P = 2.5e-06"),
            (1.25e-10, "P = 1.25e-10"),
            (0.05 / 1_000_000, "P = 5e-08"),
            (0.05 / 3, "P = 1.67e-02"),
        ],
    )
    def test_label(self, threshold, label):
        assert threshold_label(threshold) == label

    @pytest.mark.parametrize("threshold", [5e-8, 2.5e-6, 7.5e-9, 1.25e-10, 0.05 / 20])
    def test_label_reads_back_as_the_threshold(self, threshold):
        shown = float(threshold_label(threshold).removeprefix("P = "))

        assert shown == pytest.approx(threshold, rel=1e-9)


class TestThresholdLegendEntries:
    def test_one_line_swatch_per_line_in_the_lines_colour_and_style(self):
        entries = threshold_legend_entries(
            [(5e-8, "red"), (1e-5, "blue")], linestyle=":", linewidth=2.0
        )

        assert entries == [
            LegendEntry("P = 1e-05", "blue", marker="line", linestyle=":", linewidth=2),
            LegendEntry("P = 5e-08", "red", marker="line", linestyle=":", linewidth=2),
        ]

    def test_no_lines_no_entries(self):
        assert threshold_legend_entries([], linestyle="--", linewidth=1.0) == []


@pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
class TestManhattanKey:
    def test_key_is_on_by_default(self, backend_name, threshold_gwas_df):
        probe = PROBES[backend_name]

        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan(
            threshold_gwas_df
        )

        assert probe.legend_lines(fig) == [SIGNIFICANCE]
        assert probe.legend_corner(fig) == "upper right"

    def test_key_names_the_threshold_the_line_is_drawn_at(
        self, backend_name, threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan(
            threshold_gwas_df, significance_threshold=2.5e-6
        )

        assert probe.legend_lines(fig) == [("P = 2.5e-06", "#ff0000", "--")]
        assert probe.hline_levels(fig) == pytest.approx([-np.log10(2.5e-6)])

    def test_option_off_draws_the_line_without_a_key(
        self, backend_name, threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan(
            threshold_gwas_df, style=NO_KEY
        )

        assert probe.legend_lines(fig) == []
        assert len(probe.hline_levels(fig)) == 1

    def test_no_line_no_key(self, backend_name, threshold_gwas_df):
        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan(
            threshold_gwas_df, significance_threshold=None
        )

        assert PROBES[backend_name].legend_lines(fig) == []

    def test_swatch_follows_the_line_style(self, backend_name, threshold_gwas_df):
        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan(
            threshold_gwas_df, style=GenomeWideStyle(line_style=":")
        )

        assert PROBES[backend_name].legend_lines(fig) == [("P = 5e-08", "#ff0000", ":")]

    def test_categorical_panel_has_a_key(self, backend_name, threshold_gwas_df):
        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan(
            threshold_gwas_df, category_col="chr"
        )

        assert PROBES[backend_name].legend_lines(fig) == [SIGNIFICANCE]

    def test_one_entry_per_line_in_a_row(self, backend_name, threshold_gwas_df):
        probe = PROBES[backend_name]

        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan_qq(
            threshold_gwas_df, suggestive_threshold=1e-5
        )

        assert probe.legend_lines(fig, panel=0) == [
            ("P = 1e-05", "#0000ff", "--"),
            SIGNIFICANCE,
        ]
        assert probe.legend_is_row(fig, panel=0)

    def test_suggestive_line_alone_has_the_only_entry(
        self, backend_name, threshold_gwas_df
    ):
        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan_qq(
            threshold_gwas_df, significance_threshold=None, suggestive_threshold=1e-5
        )

        assert PROBES[backend_name].legend_lines(fig, panel=0) == [
            ("P = 1e-05", "#0000ff", "--")
        ]

    def test_qq_panel_has_no_key(self, backend_name, threshold_gwas_df):
        fig = ManhattanPlotter(species="human", backend=backend_name).plot_manhattan_qq(
            threshold_gwas_df
        )

        assert PROBES[backend_name].legend_lines(fig, panel=1) == []

    def test_every_stacked_panel_has_a_key(
        self, backend_name, threshold_gwas_df, second_threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = ManhattanPlotter(
            species="human", backend=backend_name
        ).plot_manhattan_stacked(
            [threshold_gwas_df, second_threshold_gwas_df], significance_threshold=1e-6
        )

        assert [probe.legend_lines(fig, panel=i) for i in range(2)] == [
            [("P = 1e-06", "#ff0000", "--")]
        ] * 2

    def test_stacked_manhattan_qq_keys_only_the_manhattan_panels(
        self, backend_name, threshold_gwas_df, second_threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = ManhattanPlotter(
            species="human", backend=backend_name
        ).plot_manhattan_qq_stacked([threshold_gwas_df, second_threshold_gwas_df])

        assert [probe.legend_lines(fig, panel=i) for i in range(4)] == [
            [SIGNIFICANCE],
            [],
            [SIGNIFICANCE],
            [],
        ]


@pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
class TestMiamiKey:
    def test_each_panel_keys_its_own_threshold(
        self, backend_name, threshold_gwas_df, second_threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = MiamiPlotter(species="human", backend=backend_name).plot_miami(
            threshold_gwas_df,
            second_threshold_gwas_df,
            top_threshold=1e-6,
            bottom_threshold=1e-4,
        )

        assert probe.legend_lines(fig, panel=0) == [("P = 1e-06", "#ff0000", "--")]
        assert probe.legend_lines(fig, panel=1) == [("P = 1e-04", "#ff0000", "--")]

    def test_key_sits_at_each_panels_outer_edge(
        self, backend_name, threshold_gwas_df, second_threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = MiamiPlotter(species="human", backend=backend_name).plot_miami(
            threshold_gwas_df, second_threshold_gwas_df
        )

        assert probe.legend_corner(fig, panel=0) == "upper right"
        assert probe.legend_corner(fig, panel=1) == "lower right"

    def test_panel_without_a_line_has_no_key(
        self, backend_name, threshold_gwas_df, second_threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = MiamiPlotter(species="human", backend=backend_name).plot_miami(
            threshold_gwas_df, second_threshold_gwas_df, top_threshold=None
        )

        assert probe.legend_lines(fig, panel=0) == []
        assert probe.legend_lines(fig, panel=1) == [SIGNIFICANCE]

    def test_option_off(
        self, backend_name, threshold_gwas_df, second_threshold_gwas_df
    ):
        probe = PROBES[backend_name]

        fig = MiamiPlotter(species="human", backend=backend_name).plot_miami(
            threshold_gwas_df, second_threshold_gwas_df, style=NO_KEY
        )

        assert [probe.legend_lines(fig, panel=i) for i in range(2)] == [[], []]

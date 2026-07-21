"""Tests for backend-neutral legend composition.

The pure builders are the deep surface: they own the bin-iteration and marker
policy that used to be duplicated across three backends. ``add_legend`` is the
single primitive each backend implements; matplotlib's is exercised directly
here, while plotly/bokeh native rendering is covered by the backend suite.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pylocuszoom.backends.composition import (
    LegendEntry,
    effect_legend_entries,
    eqtl_legend_entries,
    finemapping_legend_entries,
    ld_legend_entries,
)
from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend
from pylocuszoom.colors import (
    EFFECT_CONGRUENT_COLOR,
    EFFECT_INCONGRUENT_COLOR,
    EQTL_NEGATIVE_BINS,
    EQTL_POSITIVE_BINS,
    LD_BINS,
    LD_NA_COLOR,
    LEAD_SNP_COLOR,
    get_credible_set_color,
)


class TestLegendBuilders:
    def test_ld_entries_lead_snp_then_bins(self):
        entries = ld_legend_entries()
        assert entries[0] == LegendEntry("Lead SNP", LEAD_SNP_COLOR, marker="D")
        assert len(entries) == 1 + len(LD_BINS)
        for entry, ld_bin in zip(entries[1:], LD_BINS):
            assert entry == LegendEntry(ld_bin.label, ld_bin.color, marker="patch")

    def test_effect_entries_three_directions(self):
        entries = effect_legend_entries()
        assert [e.label for e in entries] == [
            "Same direction",
            "Opposite direction",
            "Missing effect",
        ]
        assert [e.color for e in entries] == [
            EFFECT_CONGRUENT_COLOR,
            EFFECT_INCONGRUENT_COLOR,
            LD_NA_COLOR,
        ]
        assert all(e.marker == "o" for e in entries)

    def test_eqtl_entries_split_by_effect_sign(self):
        entries = eqtl_legend_entries()
        n_pos = len(EQTL_POSITIVE_BINS)
        assert len(entries) == n_pos + len(EQTL_NEGATIVE_BINS)
        assert all(e.marker == "^" for e in entries[:n_pos])
        assert all(e.marker == "v" for e in entries[n_pos:])
        assert entries[0] == LegendEntry(
            EQTL_POSITIVE_BINS[0].label, EQTL_POSITIVE_BINS[0].color, marker="^"
        )

    def test_finemapping_entries_one_per_set(self):
        entries = finemapping_legend_entries([1, 3])
        assert [e.label for e in entries] == ["CS1", "CS3"]
        assert [e.color for e in entries] == [
            get_credible_set_color(1),
            get_credible_set_color(3),
        ]
        assert all(e.marker == "o" for e in entries)

    def test_finemapping_entries_empty(self):
        assert finemapping_legend_entries([]) == []


class TestMatplotlibAddLegend:
    def test_renders_labels_title_and_mixed_markers(self):
        fig, ax = plt.subplots()
        try:
            entries = [
                LegendEntry("Swatch", "#ff0000", marker="patch"),
                LegendEntry("Point", "#00ff00", marker="^"),
            ]
            legend = MatplotlibBackend().add_legend(ax, entries, title="r²")
            assert [t.get_text() for t in legend.get_texts()] == ["Swatch", "Point"]
            assert legend.get_title().get_text() == "r²"
        finally:
            plt.close(fig)

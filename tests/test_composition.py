"""Tests for backend-neutral legend composition.

The pure builders are the deep surface: they own the bin-iteration and marker
policy that used to be duplicated across three backends. ``add_legend`` is the
single primitive each backend implements; matplotlib's is exercised directly
here, while plotly/bokeh native rendering is covered by the backend suite.
"""

import matplotlib
import pytest

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


class _RecordingSecondaryBackend:
    """Records secondary-axis primitive calls to test overlay dispatch."""

    def __init__(self):
        self.calls = []

    def create_twin_axis(self, ax):
        self.calls.append("create_twin_axis")
        return ("secondary-handle", "y2")

    def fill_between_secondary(self, secondary, *args, **kwargs):
        self.calls.append("fill_between_secondary")

    def line_secondary(self, secondary, *args, **kwargs):
        self.calls.append("line_secondary")

    def set_secondary_ylim(self, secondary, *args, **kwargs):
        self.calls.append("set_secondary_ylim")

    def set_secondary_ylabel(self, secondary, *args, **kwargs):
        self.calls.append("set_secondary_ylabel")


class TestRecombinationOverlay:
    def test_drives_secondary_primitives_in_order(self):
        import pandas as pd

        from pylocuszoom.backends.composition import render_recombination_overlay

        backend = _RecordingSecondaryBackend()
        recomb = pd.DataFrame({"pos": [100, 200, 300], "rate": [10.0, 50.0, 20.0]})
        render_recombination_overlay(backend, "AX", recomb, start=150, end=300)
        assert backend.calls == [
            "create_twin_axis",
            "fill_between_secondary",
            "line_secondary",
            "set_secondary_ylim",
            "set_secondary_ylabel",
        ]

    def test_empty_region_skips_all_primitives(self):
        import pandas as pd

        from pylocuszoom.backends.composition import render_recombination_overlay

        backend = _RecordingSecondaryBackend()
        recomb = pd.DataFrame({"pos": [100, 200], "rate": [10.0, 50.0]})
        render_recombination_overlay(backend, "AX", recomb, start=500, end=600)
        assert backend.calls == []


class TestHeatmapHighlightCells:
    """The cell walk behind the lead-SNP outline."""

    def test_lead_snp_covers_row_left_of_diagonal_then_column_below(self):
        from pylocuszoom.backends.composition import heatmap_highlight_cells

        assert heatmap_highlight_cells(2, 5) == [
            (0, 2),
            (1, 2),
            (2, 2),
            (2, 3),
            (2, 4),
        ]

    def test_first_snp_is_diagonal_cell_plus_full_column(self):
        from pylocuszoom.backends.composition import heatmap_highlight_cells

        assert heatmap_highlight_cells(0, 3) == [(0, 0), (0, 1), (0, 2)]

    def test_last_snp_is_full_row_and_no_column(self):
        from pylocuszoom.backends.composition import heatmap_highlight_cells

        assert heatmap_highlight_cells(2, 3) == [(0, 2), (1, 2), (2, 2)]

    def test_every_cell_lies_in_the_lower_triangle(self):
        from pylocuszoom.backends.composition import heatmap_highlight_cells

        for idx in range(6):
            for x, y in heatmap_highlight_cells(idx, 6):
                assert x <= y

    @pytest.mark.parametrize(
        "snp_idx,n_snps",
        [(0, 0), (-1, 5), (5, 5), (6, 5), (0, -1)],
    )
    def test_out_of_bounds_raises(self, snp_idx, n_snps):
        from pylocuszoom.backends.composition import heatmap_highlight_cells

        with pytest.raises(ValueError, match="Invalid snp_idx"):
            heatmap_highlight_cells(snp_idx, n_snps)


class TestHeatmapHighlightRects:
    """Outline geometry in the same coordinates the heatmap was drawn in."""

    def test_cell_edges_sit_at_midpoints_between_neighbours(self):
        from pylocuszoom.backends.composition import cell_edges

        assert cell_edges([0.0, 1.0, 2.0]) == [(-0.5, 0.5), (0.5, 1.5), (1.5, 2.5)]

    def test_cell_edges_of_uneven_spacing_mirror_the_outer_gaps(self):
        from pylocuszoom.backends.composition import cell_edges

        assert cell_edges([0.0, 10.0, 12.0]) == [
            (-5.0, 5.0),
            (5.0, 11.0),
            (11.0, 13.0),
        ]

    def test_single_cell_spans_one_unit(self):
        from pylocuszoom.backends.composition import cell_edges

        assert cell_edges([7.0]) == [(6.5, 7.5)]

    def test_index_coordinates_give_unit_cells(self):
        from pylocuszoom.backends.composition import heatmap_highlight_rects

        coords = [0, 1, 2]

        assert heatmap_highlight_rects(2, coords, coords) == [
            (-0.5, 1.5, 1.0, 1.0),
            (0.5, 1.5, 1.0, 1.0),
            (1.5, 1.5, 1.0, 1.0),
        ]

    def test_genomic_x_coordinates_size_cells_by_spacing(self):
        from pylocuszoom.backends.composition import heatmap_highlight_rects

        rects = heatmap_highlight_rects(0, [1_000_000, 1_000_500], [0, 1])

        assert [(x, width) for x, _, width, _ in rects] == [
            (999_750.0, 500.0),
            (999_750.0, 500.0),
        ]

    def test_out_of_bounds_raises(self):
        from pylocuszoom.backends.composition import heatmap_highlight_rects

        with pytest.raises(ValueError, match="Invalid snp_idx"):
            heatmap_highlight_rects(3, [0, 1, 2], [0, 1, 2])

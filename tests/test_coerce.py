"""Tests for the matplotlib-vocabulary coercions in backends/_coerce.py.

Plotly and bokeh both translate out of matplotlib's vocabulary through this
module, so a change here moves both backends at once.
"""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends import _coerce


class TestPixels:
    """Figure sizes in inches become pixels at the nominal 100 dpi."""

    def test_inches_become_pixels_at_100_dpi(self):
        """An inch of figsize is 100 pixels on both axes."""
        assert _coerce.pixels((8.0, 4.0)) == (800, 400)

    def test_fractional_inches_truncate(self):
        """A fractional pixel count truncates rather than rounding."""
        assert _coerce.pixels((6.755, 3.999)) == (675, 399)


class TestNormalizeRatios:
    """Panel ratios are scaled to sum to one."""

    def test_none_passes_through(self):
        """No ratios means the caller wants the backend's even split."""
        assert _coerce.normalize_ratios(None) is None

    def test_ratios_sum_to_one(self):
        """Any relative sizes are rescaled to a unit total."""
        assert _coerce.normalize_ratios([3, 1]) == pytest.approx([0.75, 0.25])

    def test_already_normalized_ratios_are_unchanged(self):
        """Ratios that already sum to one survive the round trip."""
        assert _coerce.normalize_ratios([0.6, 0.4]) == pytest.approx([0.6, 0.4])


class TestSplitPixels:
    """A pixel extent is divided by ratio, or evenly when there are none."""

    def test_no_ratios_splits_evenly(self):
        """Without ratios every part gets the same integer share."""
        assert _coerce.split_pixels(900, None, 3) == [300, 300, 300]

    def test_even_split_truncates_the_remainder(self):
        """An extent that does not divide evenly loses the remainder."""
        assert _coerce.split_pixels(1000, None, 3) == [333, 333, 333]

    def test_ratios_drive_the_split(self):
        """Ratios divide the extent in proportion, whatever their total."""
        assert _coerce.split_pixels(800, [3, 1], 2) == [600, 200]


class TestMarkerDiameter:
    """Matplotlib marker areas become the diameters web backends expect."""

    def test_scalar_area_stays_scalar(self):
        """A single area returns a single diameter, not a list."""
        assert _coerce.marker_diameter(100.0) == 10.0

    def test_sequence_maps_element_wise(self):
        """One area per point returns one diameter per point."""
        assert _coerce.marker_diameter([100.0, 400.0]) == [10.0, 20.0]

    def test_small_markers_are_floored_at_six(self):
        """A tiny area is raised to 6 so the marker stays clickable."""
        assert _coerce.marker_diameter(1.0) == 6
        assert _coerce.marker_diameter([1.0, 400.0]) == [6, 20.0]


class TestMarkerColors:
    """One colour or one per point passes through as plain values."""

    def test_single_colour_stays_scalar(self):
        """A string colour is not expanded to one entry per point."""
        assert _coerce.marker_colors("#FF0000") == "#FF0000"

    def test_series_becomes_a_list(self):
        """A pandas Series is materialized so backends can serialize it."""
        assert _coerce.marker_colors(pd.Series(["#FF0000", "#00FF00"])) == [
            "#FF0000",
            "#00FF00",
        ]

    def test_ndarray_becomes_a_list(self):
        """An ndarray of colours is materialized the same way."""
        assert _coerce.marker_colors(np.array(["#FF0000", "#00FF00"])) == [
            "#FF0000",
            "#00FF00",
        ]


class TestPerPoint:
    """Every value is expanded to exactly one entry per point."""

    def test_scalar_expands_to_n_entries(self):
        """A single value is repeated for every point."""
        assert _coerce.per_point("#FF0000", 3) == ["#FF0000"] * 3

    def test_number_expands_to_n_entries(self):
        """A numeric scalar expands the same way a string does."""
        assert _coerce.per_point(0.5, 2) == [0.5, 0.5]

    def test_sequence_is_listed_unchanged(self):
        """A per-point sequence is listed without being resized."""
        assert _coerce.per_point(pd.Series([1, 2, 3]), 3) == [1, 2, 3]


class TestBroadcast:
    """Fill bounds are repeated across points, keeping a Series as an ndarray."""

    def test_scalar_repeats_across_points(self):
        """A single bound is repeated for every point."""
        assert _coerce.broadcast(0.0, 3) == [0.0, 0.0, 0.0]

    def test_series_stays_an_ndarray(self):
        """A Series is handed over as its ndarray, which bokeh packs as base64."""
        result = _coerce.broadcast(pd.Series([1.0, 2.0]), 2)

        assert isinstance(result, np.ndarray)
        assert list(result) == [1.0, 2.0]

    def test_plain_sequence_becomes_a_list(self):
        """A sequence that is not a Series is materialized as a list."""
        assert _coerce.broadcast((1.0, 2.0), 2) == [1.0, 2.0]

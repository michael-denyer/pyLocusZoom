"""Tests for LD color utilities."""

import math

from hypothesis import assume, given
from hypothesis import strategies as st

from pylocuszoom.colors import (
    CREDIBLE_SET_COLORS,
    EQTL_NEGATIVE_BINS,
    EQTL_POSITIVE_BINS,
    LD_BINS,
    LD_HEATMAP_COLORS,
    LD_NA_COLOR,
    LD_NA_LABEL,
    LEAD_SNP_HIGHLIGHT_COLOR,
    SECONDARY_HIGHLIGHT_COLOR,
    EQTLBin,
    LDBin,
    _find_eqtl_bin,
    get_credible_set_color,
    get_eqtl_color,
    get_ld_bin,
    get_ld_color,
    get_ld_color_palette,
)
from tests.strategies import ld_values


class TestGetLdColor:
    """Tests for get_ld_color function."""

    def test_high_ld_returns_red(self):
        """R² >= 0.8 should return red."""
        assert get_ld_color(0.8) == "#FF0000"
        assert get_ld_color(0.9) == "#FF0000"
        assert get_ld_color(1.0) == "#FF0000"

    def test_medium_high_ld_returns_orange(self):
        """0.6 <= R² < 0.8 should return orange."""
        assert get_ld_color(0.6) == "#FFA500"
        assert get_ld_color(0.7) == "#FFA500"
        assert get_ld_color(0.79) == "#FFA500"

    def test_medium_ld_returns_green(self):
        """0.4 <= R² < 0.6 should return green."""
        assert get_ld_color(0.4) == "#00CD00"
        assert get_ld_color(0.5) == "#00CD00"
        assert get_ld_color(0.59) == "#00CD00"

    def test_low_medium_ld_returns_cyan(self):
        """0.2 <= R² < 0.4 should return cyan."""
        assert get_ld_color(0.2) == "#00EEEE"
        assert get_ld_color(0.3) == "#00EEEE"
        assert get_ld_color(0.39) == "#00EEEE"

    def test_low_ld_returns_blue(self):
        """0.0 <= R² < 0.2 should return blue."""
        assert get_ld_color(0.0) == "#4169E1"
        assert get_ld_color(0.1) == "#4169E1"
        assert get_ld_color(0.19) == "#4169E1"

    def test_nan_returns_grey(self):
        """NaN values should return grey."""
        assert get_ld_color(float("nan")) == LD_NA_COLOR
        assert get_ld_color(math.nan) == LD_NA_COLOR

    def test_none_returns_grey(self):
        """None values should return grey."""
        assert get_ld_color(None) == LD_NA_COLOR

    def test_boundary_values(self):
        """Test exact boundary values."""
        # Each threshold should be in the higher bin
        assert get_ld_color(0.8) == "#FF0000"
        assert get_ld_color(0.6) == "#FFA500"
        assert get_ld_color(0.4) == "#00CD00"
        assert get_ld_color(0.2) == "#00EEEE"
        assert get_ld_color(0.0) == "#4169E1"


class TestGetLdBin:
    """Tests for get_ld_bin function."""

    def test_high_ld_bin_label(self):
        """R² >= 0.8 should return '0.8 - 1.0' label."""
        assert get_ld_bin(0.85) == "0.8 - 1.0"

    def test_medium_ld_bin_label(self):
        """0.4 <= R² < 0.6 should return '0.4 - 0.6' label."""
        assert get_ld_bin(0.5) == "0.4 - 0.6"

    def test_low_ld_bin_label(self):
        """R² < 0.2 should return '0.0 - 0.2' label."""
        assert get_ld_bin(0.1) == "0.0 - 0.2"

    def test_nan_returns_na_label(self):
        """NaN should return 'NA' label."""
        assert get_ld_bin(float("nan")) == LD_NA_LABEL

    def test_none_returns_na_label(self):
        """None should return 'NA' label."""
        assert get_ld_bin(None) == LD_NA_LABEL


class TestEqtlColors:
    """Tests for eQTL effect size color functions."""

    def test_positive_effect_high(self):
        """Positive effects >= 0.3 should return dark red."""
        assert get_eqtl_color(0.35) == "#8B1A1A"
        assert get_eqtl_color(0.4) == "#8B1A1A"  # Boundary
        assert get_eqtl_color(0.5) == "#8B1A1A"  # Above max

    def test_positive_effect_medium(self):
        """Positive effects 0.2-0.3 should return orange."""
        assert get_eqtl_color(0.25) == "#FF6600"
        assert get_eqtl_color(0.3) == "#8B1A1A"  # Boundary goes to higher bin

    def test_positive_effect_low(self):
        """Positive effects 0.1-0.2 should return light orange."""
        assert get_eqtl_color(0.15) == "#FFB347"
        assert get_eqtl_color(0.1) == "#FFB347"  # Boundary

    def test_positive_effect_below_threshold(self):
        """Small positive effects (0.0-0.1) return near-zero positive color."""
        assert get_eqtl_color(0.05) == "#FFDAB9"  # peach puff
        assert get_eqtl_color(0.0) == "#FFDAB9"

    def test_negative_effect_high(self):
        """Negative effects <= -0.3 should return dark blue."""
        assert get_eqtl_color(-0.35) == "#00008B"
        assert get_eqtl_color(-0.4) == "#00008B"  # Boundary
        assert get_eqtl_color(-0.5) == "#00008B"  # Below min

    def test_negative_effect_medium(self):
        """Negative effects -0.3 to -0.2 should return steel blue."""
        assert get_eqtl_color(-0.25) == "#4682B4"

    def test_negative_effect_low(self):
        """Negative effects -0.2 to -0.1 should return aquamarine."""
        assert get_eqtl_color(-0.15) == "#66CDAA"
        assert get_eqtl_color(-0.1) == "#66CDAA"  # Boundary

    def test_negative_effect_below_threshold(self):
        """Small negative effects (-0.1 to 0.0) return near-zero negative color."""
        assert get_eqtl_color(-0.05) == "#B0E0E6"  # powder blue

    def test_nan_returns_grey(self):
        """NaN should return grey."""
        assert get_eqtl_color(float("nan")) == LD_NA_COLOR
        assert get_eqtl_color(math.nan) == LD_NA_COLOR

    def test_none_returns_grey(self):
        """None should return grey."""
        assert get_eqtl_color(None) == LD_NA_COLOR


class TestEqtlBins:
    """Tests for eQTL bin selection."""

    def test_positive_effect_bins(self):
        """Positive effects should return correct bin labels."""
        assert _find_eqtl_bin(0.35).label == "0.3 : 0.4"
        assert _find_eqtl_bin(0.25).label == "0.2 : 0.3"
        assert _find_eqtl_bin(0.15).label == "0.1 : 0.2"

    def test_negative_effect_bins(self):
        """Negative effects should return correct bin labels."""
        assert _find_eqtl_bin(-0.35).label == "-0.4 : -0.3"
        assert _find_eqtl_bin(-0.25).label == "-0.3 : -0.2"
        assert _find_eqtl_bin(-0.15).label == "-0.2 : -0.1"

    def test_small_effects_return_correct_bin(self):
        """Near-zero effects return the appropriate near-zero bins."""
        assert _find_eqtl_bin(0.05).label == "0.0 : 0.1"
        assert _find_eqtl_bin(-0.05).label == "-0.1 : 0.0"

    def test_effects_beyond_the_boundary_land_in_the_outermost_bin(self):
        """An effect past the outermost edge takes that edge's bin."""
        assert _find_eqtl_bin(0.4).label == "0.3 : 0.4"
        assert _find_eqtl_bin(5.0).label == "0.3 : 0.4"
        assert _find_eqtl_bin(-0.4).label == "-0.4 : -0.3"
        assert _find_eqtl_bin(-5.0).label == "-0.4 : -0.3"


class TestCredibleSetColors:
    """Tests for credible set color functions."""

    def test_valid_cs_ids(self):
        """Valid cs_ids should return expected colors."""
        assert get_credible_set_color(1) == CREDIBLE_SET_COLORS[0]
        assert get_credible_set_color(2) == CREDIBLE_SET_COLORS[1]
        assert get_credible_set_color(10) == CREDIBLE_SET_COLORS[9]

    def test_cs_id_zero_returns_grey(self):
        """cs_id of 0 should return grey."""
        assert get_credible_set_color(0) == LD_NA_COLOR

    def test_cs_id_negative_returns_grey(self):
        """Negative cs_id should return grey."""
        assert get_credible_set_color(-1) == LD_NA_COLOR
        assert get_credible_set_color(-100) == LD_NA_COLOR

    def test_cs_id_cycles_after_10(self):
        """cs_id > 10 should cycle through colors."""
        assert get_credible_set_color(11) == CREDIBLE_SET_COLORS[0]
        assert get_credible_set_color(12) == CREDIBLE_SET_COLORS[1]


class TestPheWASColors:
    """Tests for PheWAS category colors."""

    def test_get_phewas_category_color(self):
        """Test PheWAS category color assignment."""
        from pylocuszoom.colors import PHEWAS_CATEGORY_COLORS, get_phewas_category_color

        # First category should return first color
        assert get_phewas_category_color(0) == PHEWAS_CATEGORY_COLORS[0]
        # Should cycle through colors
        n_colors = len(PHEWAS_CATEGORY_COLORS)
        assert get_phewas_category_color(n_colors) == PHEWAS_CATEGORY_COLORS[0]

    def test_phewas_category_palette(self):
        """Test PheWAS category color palette generation."""
        from pylocuszoom.colors import get_phewas_category_palette

        categories = ["Cardiovascular", "Metabolic", "Neurological"]
        palette = get_phewas_category_palette(categories)

        assert len(palette) == 3
        assert "Cardiovascular" in palette
        assert all(c.startswith("#") for c in palette.values())


class TestGetLdColorPalette:
    """Tests for get_ld_color_palette function."""

    def test_palette_contains_all_bins(self):
        """Palette should have all bin labels."""
        palette = get_ld_color_palette()
        for _, label, _ in LD_BINS:
            assert label in palette

    def test_palette_contains_na(self):
        """Palette should have NA label."""
        palette = get_ld_color_palette()
        assert LD_NA_LABEL in palette
        assert palette[LD_NA_LABEL] == LD_NA_COLOR

    def test_palette_colors_match(self):
        """Palette colors should match LD_BINS."""
        palette = get_ld_color_palette()
        for _, label, color in LD_BINS:
            assert palette[label] == color


# =============================================================================
# Property-Based Tests (Hypothesis)
# =============================================================================


class TestLdColorProperties:
    """Property-based tests for LD color assignment."""

    @given(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))
    def test_ld_value_always_maps_to_valid_color(self, r2):
        """Any valid R² value should return a valid hex color."""
        color = get_ld_color(r2)
        assert color.startswith("#")
        assert len(color) == 7  # #RRGGBB format

    @given(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))
    def test_ld_value_always_maps_to_valid_bin(self, r2):
        """Any valid R² value should return a known bin label."""
        bin_label = get_ld_bin(r2)
        valid_labels = [label for _, label, _ in LD_BINS] + [LD_NA_LABEL]
        assert bin_label in valid_labels

    @given(ld_values(allow_nan=True))
    def test_ld_nan_handled_gracefully(self, r2):
        """NaN values should return NA color/label."""
        if math.isnan(r2) if isinstance(r2, float) else False:
            assert get_ld_color(r2) == LD_NA_COLOR
            assert get_ld_bin(r2) == LD_NA_LABEL

    @given(
        st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
        st.floats(min_value=0.0, max_value=1.0, allow_nan=False),
    )
    def test_ld_bins_monotonic(self, r2_low, r2_high):
        """Higher R² should map to same or higher bin threshold."""
        assume(r2_low <= r2_high)
        # Get bin index (lower threshold = higher bin index in LD_BINS)
        bin_low = next(
            (i for i, (thresh, _, _) in enumerate(LD_BINS) if r2_low >= thresh),
            len(LD_BINS) - 1,
        )
        bin_high = next(
            (i for i, (thresh, _, _) in enumerate(LD_BINS) if r2_high >= thresh),
            len(LD_BINS) - 1,
        )
        # Lower index = higher threshold = higher LD bin
        assert bin_high <= bin_low

    @given(st.integers(min_value=0, max_value=100))
    def test_phewas_color_cycles(self, idx):
        """PheWAS colors should cycle without error."""
        from pylocuszoom.colors import PHEWAS_CATEGORY_COLORS, get_phewas_category_color

        color = get_phewas_category_color(idx)
        expected_idx = idx % len(PHEWAS_CATEGORY_COLORS)
        assert color == PHEWAS_CATEGORY_COLORS[expected_idx]

    @given(st.integers(min_value=1, max_value=50))
    def test_credible_set_color_cycles(self, cs_id):
        """Credible set colors should cycle without error."""
        from pylocuszoom.colors import CREDIBLE_SET_COLORS, get_credible_set_color

        color = get_credible_set_color(cs_id)
        expected_idx = (cs_id - 1) % len(CREDIBLE_SET_COLORS)
        assert color == CREDIBLE_SET_COLORS[expected_idx]

    @given(st.floats(min_value=-1.0, max_value=1.0, allow_nan=False))
    def test_eqtl_color_always_valid(self, effect):
        """Any eQTL effect should return a valid hex color."""
        color = get_eqtl_color(effect)
        assert color.startswith("#")
        assert len(color) == 7

    @given(st.floats(min_value=-1.0, max_value=1.0, allow_nan=False))
    def test_eqtl_bin_always_valid(self, effect):
        """Any eQTL effect should return a known bin label."""
        bin_label = _find_eqtl_bin(effect).label
        valid_labels = (
            [label for _, _, label, _ in EQTL_POSITIVE_BINS]
            + [label for _, _, label, _ in EQTL_NEGATIVE_BINS]
            + [LD_NA_LABEL]
        )
        assert bin_label in valid_labels

    @given(st.integers(min_value=-10, max_value=0))
    def test_invalid_cs_id_returns_na(self, cs_id):
        """Invalid cs_id (<=0) should return NA color."""
        assert get_credible_set_color(cs_id) == LD_NA_COLOR


# =============================================================================
# LD Heatmap Color Tests
# =============================================================================


class TestLdHeatmapColors:
    """Tests for LD heatmap color constants."""

    def test_ld_heatmap_colors_has_two_elements(self):
        """LD_HEATMAP_COLORS should have exactly 2 elements (start, end)."""
        assert len(LD_HEATMAP_COLORS) == 2

    def test_ld_heatmap_colors_are_valid_hex(self):
        """LD_HEATMAP_COLORS should be valid hex color codes."""
        for color in LD_HEATMAP_COLORS:
            assert color.startswith("#"), f"{color} should start with #"
            assert len(color) == 7, f"{color} should be 7 characters (#RRGGBB)"
            # Verify hex digits
            hex_part = color[1:]
            int(hex_part, 16)  # Raises ValueError if not valid hex

    def test_ld_heatmap_colors_white_to_red(self):
        """LD_HEATMAP_COLORS should be white-to-red gradient."""
        assert LD_HEATMAP_COLORS[0] == "#FFFFFF"  # white
        assert LD_HEATMAP_COLORS[1] == "#FF0000"  # red

    def test_lead_snp_highlight_color_is_valid_hex(self):
        """LEAD_SNP_HIGHLIGHT_COLOR should be a valid hex color code."""
        assert LEAD_SNP_HIGHLIGHT_COLOR.startswith("#")
        assert len(LEAD_SNP_HIGHLIGHT_COLOR) == 7
        int(LEAD_SNP_HIGHLIGHT_COLOR[1:], 16)

    def test_secondary_highlight_color_is_valid_hex(self):
        """SECONDARY_HIGHLIGHT_COLOR should be a valid hex color code."""
        assert SECONDARY_HIGHLIGHT_COLOR.startswith("#")
        assert len(SECONDARY_HIGHLIGHT_COLOR) == 7
        int(SECONDARY_HIGHLIGHT_COLOR[1:], 16)

    def test_highlight_colors_are_distinct(self):
        """LEAD_SNP_HIGHLIGHT_COLOR and SECONDARY_HIGHLIGHT_COLOR should be different."""
        assert LEAD_SNP_HIGHLIGHT_COLOR != SECONDARY_HIGHLIGHT_COLOR


# =============================================================================
# NamedTuple Tests
# =============================================================================


class TestLDBinNamedTuple:
    """Tests for LDBin NamedTuple."""

    def test_ld_bins_are_ldbin_instances(self):
        """All entries in LD_BINS should be LDBin instances."""
        for ld_bin in LD_BINS:
            assert isinstance(ld_bin, LDBin)

    def test_ldbin_named_access(self):
        """LDBin fields should be accessible by name."""
        first = LD_BINS[0]
        assert first.threshold == 0.8
        assert first.label == "0.8 - 1.0"
        assert first.color == "#FF0000"

    def test_ldbin_backward_compatible_with_tuple_indexing(self):
        """LDBin should still work with index access."""
        first = LD_BINS[0]
        assert first[0] == 0.8
        assert first[1] == "0.8 - 1.0"
        assert first[2] == "#FF0000"

    def test_ldbin_backward_compatible_with_unpacking(self):
        """LDBin should still work with tuple unpacking."""
        threshold, label, color = LD_BINS[0]
        assert threshold == 0.8
        assert label == "0.8 - 1.0"
        assert color == "#FF0000"


class TestEQTLBinNamedTuple:
    """Tests for EQTLBin NamedTuple."""

    def test_eqtl_positive_bins_are_eqtlbin_instances(self):
        """All entries in EQTL_POSITIVE_BINS should be EQTLBin instances."""
        for eqtl_bin in EQTL_POSITIVE_BINS:
            assert isinstance(eqtl_bin, EQTLBin)

    def test_eqtl_negative_bins_are_eqtlbin_instances(self):
        """All entries in EQTL_NEGATIVE_BINS should be EQTLBin instances."""
        for eqtl_bin in EQTL_NEGATIVE_BINS:
            assert isinstance(eqtl_bin, EQTLBin)

    def test_eqtlbin_named_access(self):
        """EQTLBin fields should be accessible by name."""
        first_pos = EQTL_POSITIVE_BINS[0]
        assert first_pos.min_val == 0.3
        assert first_pos.max_val == 0.4
        assert first_pos.label == "0.3 : 0.4"
        assert first_pos.color == "#8B1A1A"

    def test_eqtlbin_backward_compatible_with_unpacking(self):
        """EQTLBin should still work with tuple unpacking."""
        min_val, max_val, label, color = EQTL_POSITIVE_BINS[0]
        assert min_val == 0.3
        assert max_val == 0.4
        assert label == "0.3 : 0.4"
        assert color == "#8B1A1A"


class TestFindEqtlBinConsistency:
    """Tests for _find_eqtl_bin shared lookup."""

    def test_find_eqtl_bin_returns_eqtlbin(self):
        """_find_eqtl_bin should return an EQTLBin instance."""
        result = _find_eqtl_bin(0.35)
        assert isinstance(result, EQTLBin)

    def test_color_comes_from_the_found_bin(self):
        """get_eqtl_color returns the colour of the bin the effect lands in."""
        for effect in (0.35, 0.25, 0.15, 0.05, 0.0, -0.15, -0.25, -0.35, -0.5):
            assert get_eqtl_color(effect) == _find_eqtl_bin(effect).color

    @given(st.floats(min_value=-1.0, max_value=1.0, allow_nan=False))
    def test_color_always_comes_from_the_found_bin(self, effect):
        """For any valid effect, the colour comes from the bin it lands in."""
        assert get_eqtl_color(effect) == _find_eqtl_bin(effect).color

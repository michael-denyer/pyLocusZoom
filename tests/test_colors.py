"""Tests for LD color utilities."""

import math

from pylocuszoom.colors import (
    LD_BINS,
    LD_NA_COLOR,
    LD_NA_LABEL,
    get_ld_bin,
    get_ld_color,
    get_ld_color_palette,
)


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


class TestPheWASColors:
    """Tests for PheWAS category colors."""

    def test_get_phewas_category_color(self):
        """Test PheWAS category color assignment."""
        from pylocuszoom.colors import get_phewas_category_color, PHEWAS_CATEGORY_COLORS

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

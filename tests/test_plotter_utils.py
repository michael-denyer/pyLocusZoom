"""Tests for shared plotter utilities."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom._data import prepare_pvalue_data
from pylocuszoom._plotter_utils import add_significance_line
from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend


class TestTransformPvalues:
    """Tests for the prepare_pvalue_data utility."""

    def test_basic_transformation(self):
        """Test basic -log10 transformation."""
        df = pd.DataFrame({"p": [0.01, 0.001, 0.0001]})
        result = prepare_pvalue_data(df, "p")

        assert "neglog10p" in result.columns
        np.testing.assert_array_almost_equal(
            result["neglog10p"].values,
            [2.0, 3.0, 4.0],
            decimal=5,
        )

    def test_clipping_extreme_values(self):
        """Test that extremely small p-values are clipped to avoid -inf."""
        df = pd.DataFrame({"p": [1e-350, 1e-400]})
        result = prepare_pvalue_data(df, "p")

        # Should be clipped to 1e-300, giving -log10(1e-300) = 300
        assert np.isfinite(result["neglog10p"].iloc[0])
        assert result["neglog10p"].iloc[0] == pytest.approx(300, abs=1)

    def test_preserves_original_columns(self):
        """Test that original DataFrame columns are preserved."""
        df = pd.DataFrame({"p": [0.05], "snp": ["rs123"], "pos": [1000]})
        result = prepare_pvalue_data(df, "p")

        assert "snp" in result.columns
        assert "pos" in result.columns
        assert result["snp"].iloc[0] == "rs123"


class TestAddSignificanceLine:
    """Tests for the add_significance_line utility.

    Assertions query the rendered matplotlib axes directly per
    CLAUDE.md's observable-outputs rule.
    """

    def test_adds_line_at_threshold(self):
        """Significance line at -log10(5e-8) with red dashed style."""
        backend = MatplotlibBackend()
        fig, ax = plt.subplots()

        add_significance_line(backend, ax, 5e-8)

        lines = ax.get_lines()
        assert len(lines) == 1, "exactly one horizontal line should be drawn"
        line = lines[0]

        # Horizontal line: both y-data points equal -log10(5e-8) ≈ 7.301
        y_data = line.get_ydata()
        assert y_data[0] == pytest.approx(7.301, abs=0.01)
        assert y_data[-1] == pytest.approx(7.301, abs=0.01)

        assert line.get_color() == "red"
        assert line.get_linestyle() == "--"

    def test_skips_when_threshold_is_none(self):
        """No line is added when threshold is None."""
        backend = MatplotlibBackend()
        fig, ax = plt.subplots()

        add_significance_line(backend, ax, None)
        assert ax.get_lines() == [], "no lines should be drawn"

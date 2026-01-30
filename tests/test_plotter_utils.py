"""Tests for shared plotter utilities."""

from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

from pylocuszoom._plotter_utils import add_significance_line, transform_pvalues


class TestTransformPvalues:
    """Tests for the transform_pvalues utility."""

    def test_basic_transformation(self):
        """Test basic -log10 transformation."""
        df = pd.DataFrame({"p": [0.01, 0.001, 0.0001]})
        result = transform_pvalues(df, "p")

        assert "neglog10p" in result.columns
        np.testing.assert_array_almost_equal(
            result["neglog10p"].values,
            [2.0, 3.0, 4.0],
            decimal=5,
        )

    def test_clipping_extreme_values(self):
        """Test that extremely small p-values are clipped to avoid -inf."""
        df = pd.DataFrame({"p": [1e-350, 1e-400]})
        result = transform_pvalues(df, "p")

        # Should be clipped to 1e-300, giving -log10(1e-300) = 300
        assert np.isfinite(result["neglog10p"].iloc[0])
        assert result["neglog10p"].iloc[0] == pytest.approx(300, abs=1)

    def test_preserves_original_columns(self):
        """Test that original DataFrame columns are preserved."""
        df = pd.DataFrame({"p": [0.05], "snp": ["rs123"], "pos": [1000]})
        result = transform_pvalues(df, "p")

        assert "snp" in result.columns
        assert "pos" in result.columns
        assert result["snp"].iloc[0] == "rs123"


class TestAddSignificanceLine:
    """Tests for the add_significance_line utility."""

    def test_adds_line_at_threshold(self):
        """Test that significance line is added at correct position."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        add_significance_line(mock_backend, mock_ax, 5e-8)

        mock_backend.axhline.assert_called_once()
        call_kwargs = mock_backend.axhline.call_args[1]
        # -log10(5e-8) ≈ 7.3
        assert call_kwargs["y"] == pytest.approx(7.3, abs=0.1)
        assert call_kwargs["color"] == "red"
        assert call_kwargs["linestyle"] == "--"

    def test_skips_when_threshold_is_none(self):
        """Test that no line is added when threshold is None."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        add_significance_line(mock_backend, mock_ax, None)

        mock_backend.axhline.assert_not_called()

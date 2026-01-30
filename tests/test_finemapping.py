"""Tests for fine-mapping/SuSiE data handling."""

from unittest.mock import MagicMock

import pandas as pd
import pytest

from pylocuszoom.finemapping import (
    FinemappingValidationError,
    filter_by_credible_set,
    filter_finemapping_by_region,
    get_credible_sets,
    get_top_pip_variants,
    plot_finemapping,
    prepare_finemapping_for_plotting,
    validate_finemapping_df,
)


@pytest.fixture
def finemapping_df():
    """Create sample fine-mapping DataFrame."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "pip": [0.95, 0.02, 0.8, 0.1, 0.01],
            "cs": [1, 0, 2, 2, 0],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        }
    )


class TestValidateFinemappingDf:
    """Tests for validate_finemapping_df function."""

    def test_valid_df(self, finemapping_df):
        """Should not raise for valid DataFrame."""
        validate_finemapping_df(finemapping_df)

    def test_missing_pos_col(self, finemapping_df):
        """Should raise for missing position column."""
        df = finemapping_df.drop(columns=["pos"])
        with pytest.raises(FinemappingValidationError, match="Missing columns"):
            validate_finemapping_df(df)

    def test_missing_pip_col(self, finemapping_df):
        """Should raise for missing PIP column."""
        df = finemapping_df.drop(columns=["pip"])
        with pytest.raises(FinemappingValidationError, match="Missing columns"):
            validate_finemapping_df(df)

    def test_invalid_pip_values(self, finemapping_df):
        """Should raise for PIP values outside [0, 1]."""
        df = finemapping_df.copy()
        df.loc[0, "pip"] = 1.5
        with pytest.raises(FinemappingValidationError, match="values > 1"):
            validate_finemapping_df(df)

    def test_custom_column_names(self):
        """Should accept custom column names."""
        df = pd.DataFrame({"position": [1000, 2000], "probability": [0.5, 0.3]})
        validate_finemapping_df(df, pos_col="position", pip_col="probability")


class TestFilterFinemappingByRegion:
    """Tests for filter_finemapping_by_region function."""

    def test_filter_by_position(self, finemapping_df):
        """Should filter to region bounds."""
        result = filter_finemapping_by_region(
            finemapping_df, chrom=1, start=1500, end=4500
        )
        assert len(result) == 3
        assert set(result["pos"]) == {2000, 3000, 4000}

    def test_filter_with_chrom(self):
        """Should filter by chromosome when column exists."""
        df = pd.DataFrame(
            {
                "chr": ["1", "1", "2"],
                "pos": [1000, 2000, 1500],
                "pip": [0.5, 0.3, 0.8],
            }
        )
        result = filter_finemapping_by_region(df, chrom=1, start=0, end=3000)
        assert len(result) == 2


class TestGetCredibleSets:
    """Tests for get_credible_sets function."""

    def test_returns_unique_cs(self, finemapping_df):
        """Should return sorted unique credible set IDs."""
        result = get_credible_sets(finemapping_df)
        assert result == [1, 2]

    def test_excludes_zero(self):
        """Should exclude cs=0 (not in credible set)."""
        df = pd.DataFrame({"pos": [1, 2, 3], "pip": [0.5, 0.3, 0.2], "cs": [0, 0, 0]})
        result = get_credible_sets(df)
        assert result == []

    def test_no_cs_column(self):
        """Should return empty list if no cs column."""
        df = pd.DataFrame({"pos": [1000, 2000], "pip": [0.5, 0.3]})
        result = get_credible_sets(df)
        assert result == []


class TestFilterByCredibleSet:
    """Tests for filter_by_credible_set function."""

    def test_filter_to_cs(self, finemapping_df):
        """Should filter to specific credible set."""
        result = filter_by_credible_set(finemapping_df, cs_id=2)
        assert len(result) == 2
        assert set(result["pos"]) == {3000, 4000}

    def test_missing_cs_col(self):
        """Should raise if cs column missing."""
        df = pd.DataFrame({"pos": [1000], "pip": [0.5]})
        with pytest.raises(FinemappingValidationError):
            filter_by_credible_set(df, cs_id=1)


class TestGetTopPipVariants:
    """Tests for get_top_pip_variants function."""

    def test_returns_top_n(self, finemapping_df):
        """Should return top N by PIP."""
        result = get_top_pip_variants(finemapping_df, n=2)
        assert len(result) == 2
        assert list(result["pip"]) == [0.95, 0.8]

    def test_respects_threshold(self, finemapping_df):
        """Should filter by PIP threshold."""
        result = get_top_pip_variants(finemapping_df, n=10, pip_threshold=0.5)
        assert len(result) == 2


class TestPrepareFinemappingForPlotting:
    """Tests for prepare_finemapping_for_plotting function."""

    def test_sorts_by_position(self, finemapping_df):
        """Should sort by position."""
        # Shuffle first
        df = finemapping_df.sample(frac=1, random_state=42)
        result = prepare_finemapping_for_plotting(df)
        assert list(result["pos"]) == sorted(finemapping_df["pos"])

    def test_filters_by_region(self, finemapping_df):
        """Should filter by region when specified."""
        result = prepare_finemapping_for_plotting(
            finemapping_df, chrom=1, start=1500, end=3500
        )
        assert len(result) == 2


class TestPlotFinemapping:
    """Tests for plot_finemapping function."""

    def test_plot_finemapping_plots_pip_line(self):
        """Test that plot_finemapping renders PIP as line."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "pip": [0.1, 0.5, 0.2],
            }
        )

        plot_finemapping(mock_backend, mock_ax, df)

        # Should call line method for PIP
        mock_backend.line.assert_called_once()

    def test_plot_finemapping_colors_by_credible_set(self):
        """Test that plot_finemapping colors points by credible set."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000, 4000],
                "pip": [0.1, 0.5, 0.2, 0.05],
                "cs": [1, 1, 2, 0],
            }
        )

        plot_finemapping(mock_backend, mock_ax, df, cs_col="cs")

        # Should call scatter for each credible set
        assert mock_backend.scatter.call_count >= 2  # CS 1 and CS 2

    def test_plot_finemapping_handles_missing_cs_col(self):
        """Test that plot_finemapping works without credible set column."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "pip": [0.1, 0.5, 0.2],
            }
        )

        # Should not raise
        plot_finemapping(mock_backend, mock_ax, df, cs_col=None)

        mock_backend.line.assert_called_once()

    def test_plot_finemapping_with_pip_threshold(self):
        """Test that pip_threshold filters scatter points."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "pip": [0.05, 0.5, 0.02],  # Only 0.5 > 0.1 threshold
            }
        )

        plot_finemapping(mock_backend, mock_ax, df, cs_col=None, pip_threshold=0.1)

        # Line should still be plotted
        mock_backend.line.assert_called_once()
        # With no credible sets and pip_threshold > 0, scatter should be called for high-PIP points
        mock_backend.scatter.assert_called_once()

    def test_plot_finemapping_empty_dataframe(self):
        """Test that plot_finemapping handles empty DataFrame."""
        mock_backend = MagicMock()
        mock_ax = MagicMock()

        df = pd.DataFrame(
            {
                "pos": [],
                "pip": [],
            }
        )

        # Should not raise
        plot_finemapping(mock_backend, mock_ax, df)

        # Line should still be called (with empty data)
        mock_backend.line.assert_called_once()

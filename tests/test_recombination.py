"""Tests for recombination rate overlay module."""

from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pylocuszoom.recombination import (
    RECOMB_COLOR,
    add_recombination_overlay,
    get_default_data_dir,
    get_recombination_rate_for_region,
    load_recombination_map,
)


class TestGetDefaultDataDir:
    """Tests for get_default_data_dir function."""

    def test_returns_path_object(self):
        """Should return a Path object."""
        result = get_default_data_dir()
        assert isinstance(result, Path)

    def test_path_ends_with_recombination_maps(self):
        """Path should end with recombination_maps directory."""
        result = get_default_data_dir()
        assert result.name == "recombination_maps"

    @patch.dict("os.environ", {"DATABRICKS_RUNTIME_VERSION": "12.2"})
    @patch("os.path.exists")
    def test_uses_dbfs_on_databricks(self, mock_exists):
        """Should use /dbfs path when on Databricks."""
        mock_exists.return_value = True
        result = get_default_data_dir()
        assert "/dbfs" in str(result)


class TestAddRecombinationOverlay:
    """Tests for add_recombination_overlay function."""

    @pytest.fixture
    def sample_recomb_df(self):
        """Sample recombination rate DataFrame."""
        return pd.DataFrame(
            {
                "pos": [1000000, 1200000, 1400000, 1600000, 1800000, 2000000],
                "rate": [0.5, 1.2, 2.5, 1.8, 0.8, 0.3],
            }
        )

    def test_creates_secondary_axis(self, sample_recomb_df):
        """Should create secondary y-axis for recombination rate."""
        fig, ax = plt.subplots()
        recomb_ax = add_recombination_overlay(
            ax, sample_recomb_df, start=1000000, end=2000000
        )
        assert recomb_ax is not None
        plt.close(fig)

    def test_returns_none_for_empty_region(self):
        """Should return None when no data in region."""
        fig, ax = plt.subplots()
        empty_df = pd.DataFrame(columns=["pos", "rate"])
        recomb_ax = add_recombination_overlay(ax, empty_df, start=1000000, end=2000000)
        assert recomb_ax is None
        plt.close(fig)

    def test_filters_to_region(self, sample_recomb_df):
        """Should only plot data within the specified region."""
        fig, ax = plt.subplots()
        # Request only subset of data
        recomb_ax = add_recombination_overlay(
            ax, sample_recomb_df, start=1300000, end=1700000
        )
        assert recomb_ax is not None
        plt.close(fig)

    def test_sets_axis_label(self, sample_recomb_df):
        """Should set y-axis label for recombination rate."""
        fig, ax = plt.subplots()
        recomb_ax = add_recombination_overlay(
            ax, sample_recomb_df, start=1000000, end=2000000
        )
        ylabel = recomb_ax.get_ylabel()
        assert "Recombination" in ylabel
        plt.close(fig)

    def test_uses_correct_color(self, sample_recomb_df):
        """Should use the defined recombination color."""
        assert RECOMB_COLOR == "#7FCDFF"  # Light blue

    def test_sets_ylim_minimum(self, sample_recomb_df):
        """Y-axis should start at 0."""
        fig, ax = plt.subplots()
        recomb_ax = add_recombination_overlay(
            ax, sample_recomb_df, start=1000000, end=2000000
        )
        ylim = recomb_ax.get_ylim()
        assert ylim[0] == 0
        plt.close(fig)


class TestLoadRecombinationMap:
    """Tests for load_recombination_map function."""

    def test_raises_for_missing_file(self, tmp_path):
        """Should raise FileNotFoundError when map file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="Recombination map not found"):
            load_recombination_map(chrom=1, data_dir=str(tmp_path))

    def test_loads_valid_map_file(self, tmp_path):
        """Should load and parse valid recombination map file."""
        # Create test file
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n1\t5000\t1.2\t0.005\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        result = load_recombination_map(chrom=1, data_dir=str(tmp_path))

        assert len(result) == 2
        assert "pos" in result.columns
        assert "rate" in result.columns
        assert result["pos"].iloc[0] == 1000
        assert result["rate"].iloc[0] == 0.5

    def test_handles_chr_prefix_in_argument(self, tmp_path):
        """Should handle 'chr' prefix in chromosome argument."""
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        # Should work with "chr1" argument
        result = load_recombination_map(chrom="chr1", data_dir=str(tmp_path))
        assert len(result) == 1


class TestGetRecombinationRateForRegion:
    """Tests for get_recombination_rate_for_region function."""

    def test_filters_to_region(self, tmp_path):
        """Should return only data within specified region."""
        # Create test file with data spanning 1000-10000
        map_content = (
            "chr\tpos\trate\tcM\n"
            "1\t1000\t0.5\t0.001\n"
            "1\t3000\t1.2\t0.003\n"
            "1\t5000\t2.0\t0.005\n"
            "1\t7000\t1.5\t0.007\n"
            "1\t10000\t0.8\t0.010\n"
        )
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        result = get_recombination_rate_for_region(
            chrom=1, start=2000, end=6000, data_dir=str(tmp_path)
        )

        # Should only include positions 3000 and 5000
        assert len(result) == 2
        assert 3000 in result["pos"].values
        assert 5000 in result["pos"].values
        assert 1000 not in result["pos"].values

    def test_returns_only_pos_and_rate_columns(self, tmp_path):
        """Should return DataFrame with only pos and rate columns."""
        map_content = "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.001\n"
        map_file = tmp_path / "chr1_recomb.tsv"
        map_file.write_text(map_content)

        result = get_recombination_rate_for_region(
            chrom=1, start=0, end=2000, data_dir=str(tmp_path)
        )

        assert list(result.columns) == ["pos", "rate"]

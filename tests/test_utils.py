"""Tests for utility functions in pylocuszoom.utils."""

import pandas as pd
import pytest

from pylocuszoom.utils import filter_by_region


class TestFilterByRegion:
    """Test filter_by_region() function."""

    # Basic position filtering

    def test_basic_position_filtering(self):
        """Filter returns only rows within position bounds."""
        df = pd.DataFrame({"pos": [1000, 2000, 3000, 4000], "value": [1, 2, 3, 4]})
        result = filter_by_region(df, region=(1, 1500, 3500), pos_col="pos")

        assert len(result) == 2
        assert list(result["pos"]) == [2000, 3000]

    def test_inclusive_bounds(self):
        """Position bounds are inclusive (>= start, <= end)."""
        df = pd.DataFrame({"pos": [1000, 2000, 3000, 4000], "value": [1, 2, 3, 4]})
        result = filter_by_region(df, region=(1, 2000, 3000), pos_col="pos")

        # Both boundary values should be included
        assert len(result) == 2
        assert list(result["pos"]) == [2000, 3000]

    # Chromosome column handling

    def test_no_chromosome_column_filters_by_position_only(self):
        """When chrom_col not in DataFrame, filter by position only."""
        df = pd.DataFrame({"pos": [1000, 2000, 3000], "value": [1, 2, 3]})
        # No 'chrom' column exists
        result = filter_by_region(df, region=(1, 1500, 2500), pos_col="pos")

        # Should still filter by position
        assert len(result) == 1
        assert result["pos"].iloc[0] == 2000

    def test_chromosome_filtering_with_column_present(self):
        """When chrom_col exists, filter by both chromosome and position."""
        df = pd.DataFrame(
            {
                "chrom": [1, 1, 2, 2],
                "pos": [1000, 2000, 1000, 2000],
                "value": [1, 2, 3, 4],
            }
        )
        result = filter_by_region(df, region=(1, 500, 2500), pos_col="pos")

        # Should only return chromosome 1 rows
        assert len(result) == 2
        assert list(result["chrom"]) == [1, 1]

    # Chromosome type coercion (int vs str, chr prefix)

    def test_chromosome_type_coercion_int_to_str(self):
        """Region chrom=1 (int) matches df['chrom']='1' (str)."""
        df = pd.DataFrame({"chrom": ["1", "1", "2"], "pos": [1000, 2000, 1000]})
        result = filter_by_region(df, region=(1, 500, 2500), pos_col="pos")

        assert len(result) == 2

    def test_chromosome_type_coercion_str_to_int(self):
        """Region chrom='1' (str) matches df['chrom']=1 (int)."""
        df = pd.DataFrame({"chrom": [1, 1, 2], "pos": [1000, 2000, 1000]})
        result = filter_by_region(df, region=("1", 500, 2500), pos_col="pos")

        assert len(result) == 2

    def test_chromosome_chr_prefix_in_region(self):
        """Region chrom='chr1' matches df['chrom']='1' or df['chrom']=1."""
        df = pd.DataFrame({"chrom": [1, 1, 2], "pos": [1000, 2000, 1000]})
        result = filter_by_region(df, region=("chr1", 500, 2500), pos_col="pos")

        assert len(result) == 2

    def test_chromosome_chr_prefix_in_dataframe(self):
        """Region chrom=1 matches df['chrom']='chr1'."""
        df = pd.DataFrame(
            {"chrom": ["chr1", "chr1", "chr2"], "pos": [1000, 2000, 1000]}
        )
        result = filter_by_region(df, region=(1, 500, 2500), pos_col="pos")

        assert len(result) == 2

    def test_chromosome_x_matching(self):
        """Chromosome X matching works across type variations."""
        df = pd.DataFrame({"chrom": ["X", "X", "1"], "pos": [1000, 2000, 1000]})
        result = filter_by_region(df, region=("chrX", 500, 2500), pos_col="pos")

        assert len(result) == 2

    # Empty result handling

    def test_empty_result_region_outside_data_range(self):
        """Region outside data range returns empty DataFrame (not error)."""
        df = pd.DataFrame({"pos": [1000, 2000, 3000], "value": [1, 2, 3]})
        result = filter_by_region(df, region=(1, 5000, 6000), pos_col="pos")

        assert len(result) == 0
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ["pos", "value"]

    def test_empty_result_wrong_chromosome(self):
        """Wrong chromosome returns empty DataFrame."""
        df = pd.DataFrame({"chrom": [1, 1, 1], "pos": [1000, 2000, 3000]})
        result = filter_by_region(df, region=(2, 500, 3500), pos_col="pos")

        assert len(result) == 0
        assert isinstance(result, pd.DataFrame)

    # Returns copy, not view

    def test_returns_copy_not_view(self):
        """Modifying result does not affect original DataFrame."""
        df = pd.DataFrame({"pos": [1000, 2000, 3000], "value": [1, 2, 3]})
        result = filter_by_region(df, region=(1, 500, 2500), pos_col="pos")

        # Modify the result
        result.loc[result.index[0], "value"] = 999

        # Original should be unchanged
        assert df["value"].iloc[0] == 1
        assert df["value"].iloc[1] == 2

    # Missing position column

    def test_missing_position_column_raises_keyerror(self):
        """Missing position column raises KeyError with helpful message."""
        df = pd.DataFrame({"wrong_col": [1000, 2000], "value": [1, 2]})

        with pytest.raises(KeyError) as exc_info:
            filter_by_region(df, region=(1, 500, 2500), pos_col="pos")

        error_msg = str(exc_info.value)
        assert "pos" in error_msg
        assert "wrong_col" in error_msg or "Available" in error_msg

    # Custom column names

    def test_custom_column_names(self):
        """Custom chrom_col and pos_col parameters work."""
        df = pd.DataFrame(
            {
                "chromosome": [1, 1, 2],
                "position": [1000, 2000, 1000],
                "value": [1, 2, 3],
            }
        )
        result = filter_by_region(
            df,
            region=(1, 500, 2500),
            chrom_col="chromosome",
            pos_col="position",
        )

        assert len(result) == 2
        assert list(result["chromosome"]) == [1, 1]

"""Tests for utility functions in pylocuszoom.utils."""

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from pylocuszoom.exceptions import ValidationError
from pylocuszoom.utils import (
    filter_by_region,
    is_spark_dataframe,
    normalize_chrom,
    to_pandas,
)


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

    def test_missing_position_column_raises_validation_error(self):
        """Missing position column raises ValidationError with helpful message."""
        df = pd.DataFrame({"wrong_col": [1000, 2000], "value": [1, 2]})

        with pytest.raises(ValidationError) as exc_info:
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


class TestIsSparkDataFrame:
    """Tests for is_spark_dataframe function."""

    def test_pandas_dataframe_returns_false(self):
        """pandas DataFrame returns False."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        assert is_spark_dataframe(df) is False

    def test_dict_returns_false(self):
        """Dictionary returns False."""
        assert is_spark_dataframe({"a": 1}) is False

    def test_list_returns_false(self):
        """List returns False."""
        assert is_spark_dataframe([1, 2, 3]) is False

    def test_mock_spark_dataframe_returns_true(self):
        """Mock PySpark DataFrame returns True."""
        mock_df = MagicMock()
        mock_df.__class__.__name__ = "DataFrame"
        mock_df.__class__.__module__ = "pyspark.sql.dataframe"

        assert is_spark_dataframe(mock_df) is True


class TestToPandas:
    """Tests for to_pandas function."""

    def test_pandas_dataframe_passthrough(self):
        """pandas DataFrame is returned as-is."""
        df = pd.DataFrame({"a": [1, 2, 3]})
        result = to_pandas(df)
        assert result is df  # Same object

    def test_unsupported_type_raises_error(self):
        """Unsupported type raises TypeError."""
        with pytest.raises(TypeError, match="Unsupported DataFrame type"):
            to_pandas([1, 2, 3])

    def test_object_with_to_pandas_method(self):
        """Object with to_pandas() method uses that method."""
        expected = pd.DataFrame({"a": [1, 2, 3]})
        mock_obj = MagicMock()
        mock_obj.to_pandas.return_value = expected
        # Make sure it's not detected as Spark
        mock_obj.__class__.__name__ = "CustomDataFrame"
        mock_obj.__class__.__module__ = "custom"

        result = to_pandas(mock_obj)
        assert result.equals(expected)
        mock_obj.to_pandas.assert_called_once()

    def test_object_with_toPandas_method(self):
        """Object with toPandas() method (Spark-style) uses that method."""
        expected = pd.DataFrame({"a": [1, 2, 3]})
        mock_obj = MagicMock()
        # Remove to_pandas but have toPandas
        del mock_obj.to_pandas
        mock_obj.toPandas.return_value = expected
        mock_obj.__class__.__name__ = "CustomDataFrame"
        mock_obj.__class__.__module__ = "custom"

        result = to_pandas(mock_obj)
        assert result.equals(expected)
        mock_obj.toPandas.assert_called_once()




class TestNormalizeChrom:
    """Tests for normalize_chrom function."""

    def test_integer_input(self):
        """Integer chromosome returns string."""
        assert normalize_chrom(1) == "1"
        assert normalize_chrom(22) == "22"

    def test_string_without_prefix(self):
        """String without chr prefix returns unchanged."""
        assert normalize_chrom("1") == "1"
        assert normalize_chrom("X") == "X"

    def test_string_with_prefix(self):
        """String with chr prefix has it removed."""
        assert normalize_chrom("chr1") == "1"
        assert normalize_chrom("chrX") == "X"


class TestPlatformCacheBase:
    """One cache root for recombination maps and Ensembl annotations."""

    def test_windows_uses_localappdata(self, tmp_path, monkeypatch):
        from pylocuszoom.utils import _platform_cache_base

        monkeypatch.setattr("sys.platform", "win32")
        monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))

        assert _platform_cache_base() == tmp_path / "pylocuszoom"

    def test_windows_falls_back_to_appdata_local(self, tmp_path, monkeypatch):
        """An unset LOCALAPPDATA still resolves to the conventional location.

        Falling back to the profile root would drop a non-hidden cache directory
        there, and would not match the %LOCALAPPDATA% default it stands in for.
        """
        from pylocuszoom.utils import _platform_cache_base

        monkeypatch.setattr("sys.platform", "win32")
        monkeypatch.delenv("LOCALAPPDATA", raising=False)
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))

        assert _platform_cache_base() == tmp_path / "AppData" / "Local" / "pylocuszoom"

    def test_windows_ignores_dbfs_and_xdg(self, tmp_path, monkeypatch):
        """The Windows branch short-circuits before the Databricks and XDG checks."""
        from pylocuszoom.utils import _platform_cache_base

        monkeypatch.setattr("sys.platform", "win32")
        monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
        monkeypatch.setenv("XDG_CACHE_HOME", "/xdg")
        monkeypatch.setattr("os.path.exists", lambda p: p == "/dbfs")

        assert _platform_cache_base() == tmp_path / "pylocuszoom"

    def test_databricks_beats_xdg(self, monkeypatch):
        """/dbfs wins over XDG_CACHE_HOME on Linux."""
        from pylocuszoom.utils import _platform_cache_base

        monkeypatch.setattr("sys.platform", "linux")
        monkeypatch.setenv("XDG_CACHE_HOME", "/xdg")
        monkeypatch.setattr("os.path.exists", lambda p: p == "/dbfs")

        assert _platform_cache_base() == Path("/dbfs/FileStore/reference_data")

    def test_empty_xdg_is_treated_as_unset(self, tmp_path, monkeypatch):
        """An empty XDG_CACHE_HOME must not resolve to a relative path."""
        from pylocuszoom.utils import _platform_cache_base

        monkeypatch.setattr("sys.platform", "linux")
        monkeypatch.setattr("os.path.exists", lambda _: False)
        monkeypatch.setenv("XDG_CACHE_HOME", "")
        monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))

        assert _platform_cache_base() == tmp_path / ".cache" / "pylocuszoom"

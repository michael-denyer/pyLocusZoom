"""Tests for utility functions in pylocuszoom.utils."""

from unittest.mock import MagicMock

import pandas as pd
import pytest

from pylocuszoom.exceptions import ValidationError
from pylocuszoom.utils import (
    filter_by_region,
    is_spark_dataframe,
    normalize_chrom,
    to_pandas,
    validate_dataframe,
    validate_genes_df,
    validate_gwas_df,
    validate_plink_files,
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


class TestValidateDataframe:
    """Tests for validate_dataframe function."""

    def test_valid_dataframe_passes(self):
        """Valid DataFrame with required columns passes."""
        df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
        validate_dataframe(df, ["a", "b"], "test_df")  # Should not raise

    def test_missing_column_raises_error(self):
        """Missing required column raises ValidationError."""
        df = pd.DataFrame({"a": [1], "b": [2]})

        with pytest.raises(ValidationError, match="missing"):
            validate_dataframe(df, ["a", "c"], "test_df")

    def test_error_includes_available_columns(self):
        """Error message includes available columns."""
        df = pd.DataFrame({"x": [1], "y": [2]})

        with pytest.raises(ValidationError) as exc_info:
            validate_dataframe(df, ["z"], "test_df")

        assert "x" in str(exc_info.value) or "y" in str(exc_info.value)


class TestValidateGwasDf:
    """Tests for validate_gwas_df function."""

    def test_valid_gwas_passes(self):
        """Valid GWAS DataFrame passes."""
        df = pd.DataFrame({"ps": [1000], "p_wald": [0.01]})
        validate_gwas_df(df)  # Should not raise

    def test_custom_column_names(self):
        """Custom column names work."""
        df = pd.DataFrame({"pos": [1000], "pval": [0.01]})
        validate_gwas_df(df, pos_col="pos", p_col="pval")  # Should not raise

    def test_missing_position_raises(self):
        """Missing position column raises error."""
        df = pd.DataFrame({"p_wald": [0.01]})

        with pytest.raises(ValidationError):
            validate_gwas_df(df)

    def test_with_rs_col(self):
        """Including rs_col validates that column too."""
        df = pd.DataFrame({"ps": [1000], "p_wald": [0.01], "rs": ["rs123"]})
        validate_gwas_df(df, rs_col="rs")  # Should not raise

    def test_missing_rs_col_when_required(self):
        """Missing rs_col when specified raises error."""
        df = pd.DataFrame({"ps": [1000], "p_wald": [0.01]})

        with pytest.raises(ValidationError):
            validate_gwas_df(df, rs_col="rs")


class TestValidateGenesDf:
    """Tests for validate_genes_df function."""

    def test_valid_genes_passes(self):
        """Valid genes DataFrame passes."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                "end": [2000],
                "gene_name": ["BRCA1"],
            }
        )
        validate_genes_df(df)  # Should not raise

    def test_missing_column_raises(self):
        """Missing required column raises error."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1000],
                # Missing "end" and "gene_name"
            }
        )

        with pytest.raises(ValidationError):
            validate_genes_df(df)


class TestValidatePlinkFiles:
    """Tests for validate_plink_files function."""

    def test_valid_plink_files(self, tmp_path):
        """Valid PLINK fileset passes."""
        # Create all required files
        (tmp_path / "test.bed").touch()
        (tmp_path / "test.bim").touch()
        (tmp_path / "test.fam").touch()

        result = validate_plink_files(tmp_path / "test")
        assert result == tmp_path / "test"

    def test_missing_bed_raises(self, tmp_path):
        """Missing .bed file raises error."""
        (tmp_path / "test.bim").touch()
        (tmp_path / "test.fam").touch()

        with pytest.raises(ValidationError, match=".bed"):
            validate_plink_files(tmp_path / "test")

    def test_missing_multiple_files(self, tmp_path):
        """Missing multiple files lists all in error."""
        (tmp_path / "test.bed").touch()
        # Missing .bim and .fam

        with pytest.raises(ValidationError) as exc_info:
            validate_plink_files(tmp_path / "test")

        assert ".bim" in str(exc_info.value)
        assert ".fam" in str(exc_info.value)


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

"""Tests for the colocalization validation helper in pylocuszoom.schemas.

Exercises `coloc_plot_spec` under both dataset names to confirm required
columns, numeric types, and p-value range checks fire correctly for GWAS and
eQTL inputs.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import ValidationError
from pylocuszoom.schemas import coloc_plot_spec
from pylocuszoom.validation import check


@pytest.fixture
def valid_df():
    return pd.DataFrame(
        {
            "rs": ["rs1", "rs2", "rs3"],
            "pos": [100, 200, 300],
            "p": [1e-8, 1e-5, 1e-3],
        }
    )


class TestValidateColocGwasDf:
    def test_valid_df_passes(self, valid_df):
        check(valid_df, coloc_plot_spec("GWAS DataFrame", pos_col="pos", p_col="p"))

    def test_valid_df_without_rs_col_passes(self, valid_df):
        check(valid_df, coloc_plot_spec("GWAS DataFrame", pos_col="pos", p_col="p"))

    def test_missing_pos_col_raises(self, valid_df):
        df = valid_df.drop(columns=["pos"])
        with pytest.raises(ValidationError, match="pos"):
            check(df, coloc_plot_spec("GWAS DataFrame", pos_col="pos", p_col="p"))

    def test_missing_p_col_raises(self, valid_df):
        df = valid_df.drop(columns=["p"])
        with pytest.raises(ValidationError, match="p"):
            check(df, coloc_plot_spec("GWAS DataFrame", pos_col="pos", p_col="p"))

    def test_non_numeric_pos_raises(self, valid_df):
        df = valid_df.copy()
        df["pos"] = ["a", "b", "c"]
        with pytest.raises(ValidationError):
            check(df, coloc_plot_spec("GWAS DataFrame", pos_col="pos", p_col="p"))

    def test_p_equal_one_passes(self, valid_df):
        """p=1 is the closed upper bound of a p-value; must not raise."""
        df = valid_df.copy()
        df.loc[0, "p"] = 1.0
        check(df, coloc_plot_spec("GWAS DataFrame", pos_col="pos", p_col="p"))


class TestValidateColocEqtlDf:
    def test_valid_df_passes(self, valid_df):
        check(valid_df, coloc_plot_spec("eQTL DataFrame", pos_col="pos", p_col="p"))

    def test_error_message_identifies_eqtl(self, valid_df):
        """eQTL-specific wrapper should name the eQTL DataFrame in errors."""
        df = valid_df.drop(columns=["p"])
        with pytest.raises(ValidationError, match="eQTL"):
            check(df, coloc_plot_spec("eQTL DataFrame", pos_col="pos", p_col="p"))


class TestValidateColocDfSharedHelper:
    """Covers coloc_plot_spec paths that are independent of the name."""

    def test_df_name_propagates_to_error(self, valid_df):
        """df_name must appear in raised error messages."""
        df = valid_df.drop(columns=["pos"])
        with pytest.raises(ValidationError, match="custom-name"):
            check(df, coloc_plot_spec("custom-name", pos_col="pos", p_col="p"))

    def test_empty_df_with_required_cols_passes(self):
        """An empty DataFrame that has the right columns is legal."""
        df = pd.DataFrame(
            {"pos": pd.Series(dtype="int64"), "p": pd.Series(dtype="float64")}
        )
        check(df, coloc_plot_spec("empty", pos_col="pos", p_col="p"))

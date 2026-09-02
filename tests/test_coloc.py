"""Tests for the colocalization validation helper in pylocuszoom.schemas.

Exercises `validate_coloc_df` under both dataset names to confirm required
columns, numeric types, and p-value range checks fire correctly for GWAS and
eQTL inputs.
"""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.exceptions import ValidationError
from pylocuszoom.schemas import validate_coloc_df


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
        validate_coloc_df(
            valid_df, "GWAS DataFrame", pos_col="pos", p_col="p", rs_col="rs"
        )

    def test_valid_df_without_rs_col_passes(self, valid_df):
        validate_coloc_df(valid_df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_missing_pos_col_raises(self, valid_df):
        df = valid_df.drop(columns=["pos"])
        with pytest.raises(ValidationError, match="pos"):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_missing_p_col_raises(self, valid_df):
        df = valid_df.drop(columns=["p"])
        with pytest.raises(ValidationError, match="p"):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_missing_rs_col_raises_when_requested(self, valid_df):
        df = valid_df.drop(columns=["rs"])
        with pytest.raises(ValidationError, match="rs"):
            validate_coloc_df(
                df, "GWAS DataFrame", pos_col="pos", p_col="p", rs_col="rs"
            )

    def test_non_numeric_pos_raises(self, valid_df):
        df = valid_df.copy()
        df["pos"] = ["a", "b", "c"]
        with pytest.raises(ValidationError):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_non_numeric_p_raises(self, valid_df):
        """String p-values must raise ValidationError, not TypeError."""
        df = valid_df.copy()
        df["p"] = ["NS", "NS", "NS"]
        with pytest.raises(ValidationError):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_p_above_one_raises(self, valid_df):
        df = valid_df.copy()
        df.loc[0, "p"] = 1.5
        with pytest.raises(ValidationError, match=r"values\s*>\s*1"):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_p_equal_zero_raises(self, valid_df):
        """p=0 yields -inf after -log10; validator enforces exclusive min."""
        df = valid_df.copy()
        df.loc[0, "p"] = 0.0
        with pytest.raises(ValidationError):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_p_negative_raises(self, valid_df):
        df = valid_df.copy()
        df.loc[0, "p"] = -0.1
        with pytest.raises(ValidationError, match=r"values\s*<=\s*0"):
            validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")

    def test_p_equal_one_passes(self, valid_df):
        """p=1 is the closed upper bound of a p-value; must not raise."""
        df = valid_df.copy()
        df.loc[0, "p"] = 1.0
        validate_coloc_df(df, "GWAS DataFrame", pos_col="pos", p_col="p")


class TestValidateColocEqtlDf:
    def test_valid_df_passes(self, valid_df):
        validate_coloc_df(
            valid_df, "eQTL DataFrame", pos_col="pos", p_col="p", rs_col="rs"
        )

    def test_error_message_identifies_eqtl(self, valid_df):
        """eQTL-specific wrapper should name the eQTL DataFrame in errors."""
        df = valid_df.drop(columns=["p"])
        with pytest.raises(ValidationError, match="eQTL"):
            validate_coloc_df(df, "eQTL DataFrame", pos_col="pos", p_col="p")


class TestValidateColocDfSharedHelper:
    """Covers validate_coloc_df paths that are independent of the name."""

    def test_df_name_propagates_to_error(self, valid_df):
        """df_name must appear in raised error messages."""
        df = valid_df.drop(columns=["pos"])
        with pytest.raises(ValidationError, match="custom-name"):
            validate_coloc_df(df, "custom-name", pos_col="pos", p_col="p")

    def test_empty_df_with_required_cols_passes(self):
        """An empty DataFrame that has the right columns is legal."""
        df = pd.DataFrame(
            {"pos": pd.Series(dtype="int64"), "p": pd.Series(dtype="float64")}
        )
        validate_coloc_df(df, "empty", pos_col="pos", p_col="p")

    def test_nan_in_p_rejected(self, valid_df):
        """NaN p-values must be rejected, not silently dropped."""
        df = valid_df.copy()
        df.loc[0, "p"] = np.nan
        with pytest.raises(ValidationError):
            validate_coloc_df(df, "gwas", pos_col="pos", p_col="p")

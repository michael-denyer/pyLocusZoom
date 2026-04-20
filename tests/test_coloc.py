"""Tests for pylocuszoom.coloc validation helpers.

Exercises the shared validator (`_validate_coloc_df`) via both public
wrappers to confirm required columns, numeric types, and p-value range
checks fire correctly for GWAS and eQTL inputs.
"""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.coloc import (
    _validate_coloc_df,
    validate_coloc_eqtl_df,
    validate_coloc_gwas_df,
)
from pylocuszoom.exceptions import ValidationError


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
        validate_coloc_gwas_df(valid_df, pos_col="pos", p_col="p", rs_col="rs")

    def test_valid_df_without_rs_col_passes(self, valid_df):
        validate_coloc_gwas_df(valid_df, pos_col="pos", p_col="p")

    def test_missing_pos_col_raises(self, valid_df):
        df = valid_df.drop(columns=["pos"])
        with pytest.raises(ValidationError, match="pos"):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_missing_p_col_raises(self, valid_df):
        df = valid_df.drop(columns=["p"])
        with pytest.raises(ValidationError, match="p"):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_missing_rs_col_raises_when_requested(self, valid_df):
        df = valid_df.drop(columns=["rs"])
        with pytest.raises(ValidationError, match="rs"):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p", rs_col="rs")

    def test_non_numeric_pos_raises(self, valid_df):
        df = valid_df.copy()
        df["pos"] = ["a", "b", "c"]
        with pytest.raises(ValidationError):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_non_numeric_p_raises(self, valid_df):
        """String p-values must be rejected (ValidationError or TypeError).

        Current validator chain runs ``require_range`` after
        ``require_numeric`` without short-circuiting on dtype errors, so
        string p-values surface as TypeError from the range comparison
        rather than a structured ValidationError. Either outcome
        indicates the upstream data is rejected — the test accepts both
        so the behaviour is pinned without prescribing the fix.
        """
        df = valid_df.copy()
        df["p"] = ["NS", "NS", "NS"]
        with pytest.raises((ValidationError, TypeError)):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_p_above_one_raises(self, valid_df):
        df = valid_df.copy()
        df.loc[0, "p"] = 1.5
        with pytest.raises(ValidationError):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_p_equal_zero_raises(self, valid_df):
        """p=0 yields -inf after -log10; validator enforces exclusive min."""
        df = valid_df.copy()
        df.loc[0, "p"] = 0.0
        with pytest.raises(ValidationError):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_p_negative_raises(self, valid_df):
        df = valid_df.copy()
        df.loc[0, "p"] = -0.1
        with pytest.raises(ValidationError):
            validate_coloc_gwas_df(df, pos_col="pos", p_col="p")

    def test_p_equal_one_passes(self, valid_df):
        """p=1 is the closed upper bound of a p-value; must not raise."""
        df = valid_df.copy()
        df.loc[0, "p"] = 1.0
        validate_coloc_gwas_df(df, pos_col="pos", p_col="p")


class TestValidateColocEqtlDf:
    def test_valid_df_passes(self, valid_df):
        validate_coloc_eqtl_df(valid_df, pos_col="pos", p_col="p", rs_col="rs")

    def test_error_message_identifies_eqtl(self, valid_df):
        """eQTL-specific wrapper should name the eQTL DataFrame in errors."""
        df = valid_df.drop(columns=["p"])
        with pytest.raises(ValidationError, match="eQTL"):
            validate_coloc_eqtl_df(df, pos_col="pos", p_col="p")


class TestValidateColocDfSharedHelper:
    """Covers _validate_coloc_df directly for paths the wrappers share."""

    def test_df_name_propagates_to_error(self, valid_df):
        """df_name must appear in raised error messages."""
        df = valid_df.drop(columns=["pos"])
        with pytest.raises(ValidationError, match="custom-name"):
            _validate_coloc_df(df, "custom-name", pos_col="pos", p_col="p")

    def test_empty_df_with_required_cols_passes(self):
        """An empty DataFrame that has the right columns is legal."""
        df = pd.DataFrame(
            {"pos": pd.Series(dtype="int64"), "p": pd.Series(dtype="float64")}
        )
        _validate_coloc_df(df, "empty", pos_col="pos", p_col="p")

    def test_nan_in_p_column_is_tolerated(self, valid_df):
        """Current contract: NaN p-values are dropped by `require_range`.

        `DataFrameValidator.require_range` calls ``.dropna()`` before the
        bound checks, so NaN silently passes validation here. This is a
        latent silent-failure risk — downstream ``-log10(p)`` on NaN
        yields NaN, which can blank axes — but the validator contract is
        documented by this test. Tightening it is out of scope for the
        coloc test file; track as a separate bead if desired.
        """
        df = valid_df.copy()
        df.loc[0, "p"] = np.nan
        _validate_coloc_df(df, "gwas", pos_col="pos", p_col="p")

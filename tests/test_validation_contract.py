"""Prose-independent contract for the loader schema validators.

Pins what callers can rely on across a change of validation implementation:
which inputs are accepted, which raise, the exception class, and which column
names the error names. Deliberately asserts nothing about error phrasing, so
the same file passes before and after ``schemas.py`` moves onto
``DataFrameValidator``.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import LoaderValidationError
from pylocuszoom.schemas import (
    validate_eqtl_dataframe,
    validate_finemapping_dataframe,
    validate_genes_dataframe,
    validate_gwas_dataframe,
)


def _genes(**overrides):
    base = {"chr": ["1"], "start": [1000], "end": [2000], "gene_name": ["GENE1"]}
    return pd.DataFrame({**base, **overrides})


ACCEPTED = [
    (
        "gwas-minimal",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [100], "p_wald": [0.5]}),
    ),
    (
        "gwas-p-exactly-one",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [1], "p_wald": [1.0]}),
    ),
    (
        "gwas-tiny-p",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [42], "p_wald": [1e-300]}),
    ),
    (
        "gwas-extra-columns",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [100], "p_wald": [0.5], "beta": [0.3], "rs": ["rs1"]}),
    ),
    (
        "eqtl-minimal",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [100], "p_value": [0.05], "gene": ["BRCA1"]}),
    ),
    (
        "eqtl-p-exactly-one",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [100], "p_value": [1.0], "gene": ["BRCA1"]}),
    ),
    (
        "eqtl-null-p-tolerated",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [100, 200], "p_value": [0.05, None], "gene": ["A", "B"]}),
    ),
    (
        "finemapping-pip-bounds-inclusive",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": [100, 200], "pip": [0.0, 1.0]}),
    ),
    (
        "finemapping-null-pip-tolerated",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": [100, 200], "pip": [0.5, None]}),
    ),
    ("genes-minimal", validate_genes_dataframe, _genes()),
    ("genes-start-at-zero", validate_genes_dataframe, _genes(start=[0])),
    ("genes-point-feature", validate_genes_dataframe, _genes(start=[1000], end=[1000])),
]


@pytest.mark.parametrize(
    ("validator", "df"),
    [pytest.param(v, df, id=name) for name, v, df in ACCEPTED],
)
def test_accepted_input_returns_the_same_object(validator, df):
    """Accepted frames pass through by identity, never a copy."""
    assert validator(df) is df


REJECTED = [
    (
        "gwas-missing-pos",
        validate_gwas_dataframe,
        pd.DataFrame({"p_wald": [0.5]}),
        ["ps"],
    ),
    (
        "gwas-missing-both",
        validate_gwas_dataframe,
        pd.DataFrame({"chr": [1]}),
        ["ps", "p_wald"],
    ),
    (
        "gwas-position-zero",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [0], "p_wald": [0.5]}),
        ["ps"],
    ),
    (
        "gwas-position-negative",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [-1], "p_wald": [0.5]}),
        ["ps"],
    ),
    (
        "gwas-null-position",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [100, None], "p_wald": [0.5, 0.5]}),
        ["ps"],
    ),
    (
        "gwas-null-pvalue",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [100, 200], "p_wald": [0.05, None]}),
        ["p_wald"],
    ),
    (
        "gwas-pvalue-zero",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [100], "p_wald": [0.0]}),
        ["p_wald"],
    ),
    (
        "gwas-pvalue-above-one",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": [100], "p_wald": [1.5]}),
        ["p_wald"],
    ),
    (
        "gwas-both-non-numeric",
        validate_gwas_dataframe,
        pd.DataFrame({"ps": ["chr1:100"], "p_wald": ["not_a_number"]}),
        ["ps", "p_wald"],
    ),
    (
        "eqtl-missing-all",
        validate_eqtl_dataframe,
        pd.DataFrame({"extra": [1]}),
        ["pos", "p_value", "gene"],
    ),
    (
        "eqtl-missing-gene-only",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [100], "p_value": [0.05]}),
        ["gene"],
    ),
    (
        "eqtl-non-numeric-pos",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": ["chr1:100"], "p_value": [0.05], "gene": ["A"]}),
        ["pos"],
    ),
    (
        "eqtl-non-numeric-pvalue",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [100], "p_value": ["significant"], "gene": ["A"]}),
        ["p_value"],
    ),
    (
        "eqtl-pvalue-zero",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [100], "p_value": [0.0], "gene": ["A"]}),
        ["p_value"],
    ),
    (
        "eqtl-non-positive-pos",
        validate_eqtl_dataframe,
        pd.DataFrame({"pos": [-5, 100], "p_value": [0.05, 0.1], "gene": ["A", "B"]}),
        ["pos"],
    ),
    (
        "finemapping-missing-both",
        validate_finemapping_dataframe,
        pd.DataFrame({"other": [1]}),
        ["pos", "pip"],
    ),
    (
        "finemapping-non-numeric-pip",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": [100], "pip": ["high"]}),
        ["pip"],
    ),
    (
        "finemapping-non-numeric-pos",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": ["abc"], "pip": [0.5]}),
        ["pos"],
    ),
    (
        "finemapping-pip-above-one",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": [100], "pip": [1.5]}),
        ["pip"],
    ),
    (
        "finemapping-pip-negative",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": [100], "pip": [-0.1]}),
        ["pip"],
    ),
    (
        "finemapping-non-positive-pos",
        validate_finemapping_dataframe,
        pd.DataFrame({"pos": [0], "pip": [0.5]}),
        ["pos"],
    ),
    (
        "genes-missing-all",
        validate_genes_dataframe,
        pd.DataFrame({"extra": [1]}),
        ["chr", "start", "end", "gene_name"],
    ),
    (
        "genes-non-numeric-start",
        validate_genes_dataframe,
        _genes(start=["abc"]),
        ["start"],
    ),
    ("genes-non-numeric-end", validate_genes_dataframe, _genes(end=["xyz"]), ["end"]),
    ("genes-negative-start", validate_genes_dataframe, _genes(start=[-1]), ["start"]),
    (
        "genes-end-before-start",
        validate_genes_dataframe,
        _genes(start=[2000], end=[1000]),
        ["start", "end"],
    ),
]


@pytest.mark.parametrize(
    ("validator", "df", "named_columns"),
    [pytest.param(v, df, cols, id=name) for name, v, df, cols in REJECTED],
)
def test_rejected_input_raises_and_names_the_offending_columns(
    validator, df, named_columns
):
    """Rejection is a LoaderValidationError whose message names each bad column."""
    with pytest.raises(LoaderValidationError) as exc_info:
        validator(df)

    message = str(exc_info.value)
    for column in named_columns:
        assert column in message, f"{column!r} missing from: {message}"


def test_custom_gwas_column_names_are_honoured():
    """pos_col and p_col redirect validation to the named columns."""
    df = pd.DataFrame({"position": [100], "pvalue": [0.5]})
    assert validate_gwas_dataframe(df, pos_col="position", p_col="pvalue") is df

    with pytest.raises(LoaderValidationError) as exc_info:
        validate_gwas_dataframe(df, pos_col="missing_pos", p_col="pvalue")
    assert "missing_pos" in str(exc_info.value)


def test_all_problems_are_reported_in_one_pass():
    """Validation accumulates rather than raising on the first fault."""
    df = pd.DataFrame({"ps": [0, 100], "p_wald": [0.5, 2.0]})
    with pytest.raises(LoaderValidationError) as exc_info:
        validate_gwas_dataframe(df)

    message = str(exc_info.value)
    assert "ps" in message
    assert "p_wald" in message


class TestPValueDomainHasOneOwner:
    """Strict validation and plot-time intake share the p-value domain."""

    @staticmethod
    def _accepted_by_strict(value):
        from pylocuszoom.exceptions import LoaderValidationError

        df = pd.DataFrame({"ps": [100], "p_wald": [value]})
        try:
            validate_gwas_dataframe(df)
        except LoaderValidationError:
            return False
        return True

    @staticmethod
    def _kept_by_intake(value, *, allow_zero):
        from pylocuszoom._data import prepare_pvalue_data

        df = pd.DataFrame({"p": [value]})
        return len(prepare_pvalue_data(df, "p", allow_zero=allow_zero)) == 1

    @pytest.mark.parametrize(
        "value", [1e-300, 1e-8, 0.05, 0.5, 0.999, 1.0, 1.0000001, 1.5, 2.0, -0.1]
    )
    def test_strict_agrees_with_intake_when_zero_is_disallowed(self, value):
        """Away from zero, both ends of the library accept exactly the same values."""
        assert self._accepted_by_strict(value) == self._kept_by_intake(
            value, allow_zero=False
        )

    def test_zero_is_the_only_deliberate_divergence(self):
        """p == 0 is rejected at load and kept by the Manhattan convention."""
        assert self._accepted_by_strict(0.0) is False
        assert self._kept_by_intake(0.0, allow_zero=False) is False
        assert self._kept_by_intake(0.0, allow_zero=True) is True

    def test_upper_bound_has_a_single_source(self):
        """require_pvalue and prepare_pvalue_data read the same constant."""
        from pylocuszoom._data import P_VALUE_MAX
        from pylocuszoom.exceptions import LoaderValidationError

        just_over = pd.DataFrame({"ps": [100], "p_wald": [P_VALUE_MAX * 1.000001]})
        with pytest.raises(LoaderValidationError):
            validate_gwas_dataframe(just_over)
        assert self._accepted_by_strict(P_VALUE_MAX) is True

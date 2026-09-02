"""Prose-independent contract for the loader schema validators.

Pins what callers can rely on: which inputs are accepted, which raise, the
exception class, and which column names the error names. Deliberately asserts
nothing about error phrasing, so a change of validation implementation is a
checkable claim rather than a promise. It held unchanged across the move onto
the declarative ``ColumnSpec`` engine.

This is the strict load-time tier. Plot-time intake is deliberately more
permissive; see the two-tier note in ``CONTEXT.md``.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import LoaderValidationError
from pylocuszoom.schemas import Family, Tier, spec
from pylocuszoom.validation import check


def _genes(**overrides):
    base = {"chr": ["1"], "start": [1000], "end": [2000], "gene_name": ["GENE1"]}
    return pd.DataFrame({**base, **overrides})


ACCEPTED = [
    (
        "gwas-minimal",
        Family.GWAS,
        pd.DataFrame({"pos": [100], "p_value": [0.5]}),
    ),
    (
        "gwas-p-exactly-one",
        Family.GWAS,
        pd.DataFrame({"pos": [1], "p_value": [1.0]}),
    ),
    (
        "gwas-tiny-p",
        Family.GWAS,
        pd.DataFrame({"pos": [42], "p_value": [1e-300]}),
    ),
    (
        "gwas-extra-columns",
        Family.GWAS,
        pd.DataFrame({"pos": [100], "p_value": [0.5], "beta": [0.3], "rs": ["rs1"]}),
    ),
    (
        "eqtl-minimal",
        Family.EQTL,
        pd.DataFrame({"pos": [100], "p_value": [0.05], "gene": ["BRCA1"]}),
    ),
    (
        "eqtl-p-exactly-one",
        Family.EQTL,
        pd.DataFrame({"pos": [100], "p_value": [1.0], "gene": ["BRCA1"]}),
    ),
    (
        "finemapping-pip-bounds-inclusive",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [100, 200], "pip": [0.0, 1.0]}),
    ),
    ("genes-minimal", Family.GENES, _genes()),
    ("genes-start-at-zero", Family.GENES, _genes(start=[0])),
    ("genes-point-feature", Family.GENES, _genes(start=[1000], end=[1000])),
]


@pytest.mark.parametrize(
    ("family", "df"),
    [pytest.param(f, df, id=name) for name, f, df in ACCEPTED],
)
def test_accepted_input_raises_nothing(family, df):
    """Accepted frames satisfy the load-time contract."""
    check(df, spec(family, Tier.LOAD))


REJECTED = [
    (
        "gwas-missing-pos",
        Family.GWAS,
        pd.DataFrame({"p_value": [0.5]}),
        ["pos"],
    ),
    (
        "gwas-missing-both",
        Family.GWAS,
        pd.DataFrame({"chr": [1]}),
        ["pos", "p_value"],
    ),
    (
        "gwas-position-zero",
        Family.GWAS,
        pd.DataFrame({"pos": [0], "p_value": [0.5]}),
        ["pos"],
    ),
    (
        "gwas-position-negative",
        Family.GWAS,
        pd.DataFrame({"pos": [-1], "p_value": [0.5]}),
        ["pos"],
    ),
    (
        "gwas-null-position",
        Family.GWAS,
        pd.DataFrame({"pos": [100, None], "p_value": [0.5, 0.5]}),
        ["pos"],
    ),
    (
        "gwas-null-pvalue",
        Family.GWAS,
        pd.DataFrame({"pos": [100, 200], "p_value": [0.05, None]}),
        ["p_value"],
    ),
    (
        "gwas-pvalue-zero",
        Family.GWAS,
        pd.DataFrame({"pos": [100], "p_value": [0.0]}),
        ["p_value"],
    ),
    (
        "gwas-pvalue-above-one",
        Family.GWAS,
        pd.DataFrame({"pos": [100], "p_value": [1.5]}),
        ["p_value"],
    ),
    (
        "gwas-both-non-numeric",
        Family.GWAS,
        pd.DataFrame({"pos": ["chr1:100"], "p_value": ["not_a_number"]}),
        ["pos", "p_value"],
    ),
    (
        "eqtl-missing-all",
        Family.EQTL,
        pd.DataFrame({"extra": [1]}),
        ["pos", "p_value", "gene"],
    ),
    (
        "eqtl-missing-gene-only",
        Family.EQTL,
        pd.DataFrame({"pos": [100], "p_value": [0.05]}),
        ["gene"],
    ),
    (
        "eqtl-non-numeric-pos",
        Family.EQTL,
        pd.DataFrame({"pos": ["chr1:100"], "p_value": [0.05], "gene": ["A"]}),
        ["pos"],
    ),
    (
        "eqtl-non-numeric-pvalue",
        Family.EQTL,
        pd.DataFrame({"pos": [100], "p_value": ["significant"], "gene": ["A"]}),
        ["p_value"],
    ),
    (
        "eqtl-pvalue-zero",
        Family.EQTL,
        pd.DataFrame({"pos": [100], "p_value": [0.0], "gene": ["A"]}),
        ["p_value"],
    ),
    (
        "eqtl-non-positive-pos",
        Family.EQTL,
        pd.DataFrame({"pos": [-5, 100], "p_value": [0.05, 0.1], "gene": ["A", "B"]}),
        ["pos"],
    ),
    (
        "finemapping-missing-both",
        Family.FINEMAPPING,
        pd.DataFrame({"other": [1]}),
        ["pos", "pip"],
    ),
    (
        "finemapping-non-numeric-pip",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [100], "pip": ["high"]}),
        ["pip"],
    ),
    (
        "finemapping-non-numeric-pos",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": ["abc"], "pip": [0.5]}),
        ["pos"],
    ),
    (
        "finemapping-pip-above-one",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [100], "pip": [1.5]}),
        ["pip"],
    ),
    (
        "finemapping-pip-negative",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [100], "pip": [-0.1]}),
        ["pip"],
    ),
    (
        "finemapping-non-positive-pos",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [0], "pip": [0.5]}),
        ["pos"],
    ),
    (
        "eqtl-null-pvalue",
        Family.EQTL,
        pd.DataFrame({"pos": [100, 200], "p_value": [0.05, None], "gene": ["A", "B"]}),
        ["p_value"],
    ),
    (
        "eqtl-null-pos",
        Family.EQTL,
        pd.DataFrame({"pos": [100, None], "p_value": [0.05, 0.1], "gene": ["A", "B"]}),
        ["pos"],
    ),
    (
        "finemapping-null-pip",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [100, 200], "pip": [0.5, None]}),
        ["pip"],
    ),
    (
        "finemapping-null-pos",
        Family.FINEMAPPING,
        pd.DataFrame({"pos": [100, None], "pip": [0.5, 0.6]}),
        ["pos"],
    ),
    (
        "genes-null-start",
        Family.GENES,
        _genes(
            chr=["1", "1"], start=[1000, None], end=[2000, 3000], gene_name=["A", "B"]
        ),
        ["start"],
    ),
    (
        "genes-missing-all",
        Family.GENES,
        pd.DataFrame({"extra": [1]}),
        ["chr", "start", "end", "gene_name"],
    ),
    (
        "genes-non-numeric-start",
        Family.GENES,
        _genes(start=["abc"]),
        ["start"],
    ),
    ("genes-non-numeric-end", Family.GENES, _genes(end=["xyz"]), ["end"]),
    ("genes-negative-start", Family.GENES, _genes(start=[-1]), ["start"]),
    (
        "genes-end-before-start",
        Family.GENES,
        _genes(start=[2000], end=[1000]),
        ["start", "end"],
    ),
]


@pytest.mark.parametrize(
    ("family", "df", "named_columns"),
    [pytest.param(f, df, cols, id=name) for name, f, df, cols in REJECTED],
)
def test_rejected_input_raises_and_names_the_offending_columns(
    family, df, named_columns
):
    """Rejection is a LoaderValidationError whose message names each bad column."""
    with pytest.raises(LoaderValidationError) as exc_info:
        check(df, spec(family, Tier.LOAD))

    message = str(exc_info.value)
    for column in named_columns:
        assert column in message, f"{column!r} missing from: {message}"


def test_custom_gwas_column_names_are_honoured():
    """pos_col and p_col redirect validation to the named columns."""
    df = pd.DataFrame({"position": [100], "pvalue": [0.5]})
    check(df, spec(Family.GWAS, Tier.LOAD, pos_col="position", p_col="pvalue"))

    with pytest.raises(LoaderValidationError) as exc_info:
        check(df, spec(Family.GWAS, Tier.LOAD, pos_col="missing_pos", p_col="pvalue"))
    assert "missing_pos" in str(exc_info.value)


def test_all_problems_are_reported_in_one_pass():
    """Validation accumulates rather than raising on the first fault."""
    df = pd.DataFrame({"pos": [0, 100], "p_value": [0.5, 2.0]})
    with pytest.raises(LoaderValidationError) as exc_info:
        check(df, spec(Family.GWAS, Tier.LOAD))

    message = str(exc_info.value)
    assert "pos" in message
    assert "p_value" in message


class TestPValueDomainHasOneOwner:
    """Strict validation and plot-time intake share the p-value domain."""

    @staticmethod
    def _accepted_by_strict(value):
        from pylocuszoom.exceptions import LoaderValidationError

        df = pd.DataFrame({"pos": [100], "p_value": [value]})
        try:
            check(df, spec(Family.GWAS, Tier.LOAD))
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
        """The p-value spec rule and prepare_pvalue_data read the same constant."""
        from pylocuszoom._data import P_VALUE_MAX
        from pylocuszoom.exceptions import LoaderValidationError

        just_over = pd.DataFrame({"pos": [100], "p_value": [P_VALUE_MAX * 1.000001]})
        with pytest.raises(LoaderValidationError):
            check(just_over, spec(Family.GWAS, Tier.LOAD))
        assert self._accepted_by_strict(P_VALUE_MAX) is True

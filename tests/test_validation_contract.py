"""The load-time contract for every loader schema, as one table per verdict.

Pins what callers can rely on: which inputs are accepted, which raise, the
exception class, and which column names the error names. It is the only home
for these rules; adding a column rule to ``schemas.py`` means adding a row
here.

Error phrasing is pinned only where a caller reads it. A ``REJECTED`` row
carries the columns the message must name, and optionally one fragment of the
message itself. Name the columns when any wording will do, so a change of
validation implementation stays a checkable claim rather than a promise. Add
the fragment when the count or the comparison in the message is what tells a
user which of their rows is wrong.

This is the strict load-time tier. Plot-time intake is deliberately more
permissive; see the two-tier note in ``CONTEXT.md``.
"""

import pandas as pd
import pytest

from pylocuszoom.exceptions import LoaderValidationError
from pylocuszoom.schemas import EQTL_LOAD, FINEMAPPING_LOAD, GENES_LOAD, gwas_load_spec
from pylocuszoom.validation import check


def _genes(**overrides):
    base = {"chr": ["1"], "start": [1000], "end": [2000], "gene_name": ["GENE1"]}
    return pd.DataFrame({**base, **overrides})


ACCEPTED = [
    (
        "gwas-minimal",
        gwas_load_spec(),
        pd.DataFrame({"pos": [100], "p_value": [0.5]}),
    ),
    (
        "gwas-p-exactly-one",
        gwas_load_spec(),
        pd.DataFrame({"pos": [1], "p_value": [1.0]}),
    ),
    (
        "gwas-tiny-p",
        gwas_load_spec(),
        pd.DataFrame({"pos": [42], "p_value": [1e-300]}),
    ),
    (
        "gwas-extra-columns",
        gwas_load_spec(),
        pd.DataFrame({"pos": [100], "p_value": [0.5], "beta": [0.3], "rs": ["rs1"]}),
    ),
    (
        "eqtl-minimal",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100], "p_value": [0.05], "gene": ["BRCA1"]}),
    ),
    (
        "eqtl-p-exactly-one",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100], "p_value": [1.0], "gene": ["BRCA1"]}),
    ),
    (
        "finemapping-pip-bounds-inclusive",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [100, 200], "pip": [0.0, 1.0]}),
    ),
    (
        "gwas-multi-row",
        gwas_load_spec(),
        pd.DataFrame(
            {
                "pos": [1000000, 1001000, 1002000],
                "p_value": [0.01, 0.001, 1e-8],
                "rs": ["rs1", "rs2", "rs3"],
            }
        ),
    ),
    (
        "eqtl-with-effect-column",
        EQTL_LOAD,
        pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
                "gene": ["BRCA1", "BRCA1"],
                "effect": [0.5, -0.3],
            }
        ),
    ),
    (
        "finemapping-with-credible-sets",
        FINEMAPPING_LOAD,
        pd.DataFrame(
            {
                "pos": [1000000, 1001000, 1002000],
                "pip": [0.85, 0.12, 0.03],
                "cs": [1, 1, 0],
            }
        ),
    ),
    ("genes-minimal", GENES_LOAD, _genes()),
    ("genes-start-at-zero", GENES_LOAD, _genes(start=[0])),
    ("genes-point-feature", GENES_LOAD, _genes(start=[1000], end=[1000])),
    (
        "genes-multi-row",
        GENES_LOAD,
        _genes(
            chr=["1", "1", "1"],
            start=[1000000, 1050000, 1100000],
            end=[1020000, 1080000, 1150000],
            gene_name=["GENE1", "GENE2", "GENE3"],
        ),
    ),
]


@pytest.mark.parametrize(
    ("contract", "df"),
    [pytest.param(f, df, id=name) for name, f, df in ACCEPTED],
)
def test_accepted_input_raises_nothing(contract, df):
    """Accepted frames satisfy the load-time contract."""
    check(df, contract)


REJECTED = [
    (
        "gwas-missing-pos",
        gwas_load_spec(),
        pd.DataFrame({"p_value": [0.5]}),
        ["pos"],
        "Missing columns",
    ),
    (
        "gwas-missing-pvalue",
        gwas_load_spec(),
        pd.DataFrame({"pos": [1000000, 1001000], "rs": ["rs1", "rs2"]}),
        ["p_value"],
        "Missing columns",
    ),
    (
        "gwas-missing-both",
        gwas_load_spec(),
        pd.DataFrame({"chr": [1]}),
        ["pos", "p_value"],
        "Missing columns",
    ),
    (
        "gwas-position-zero",
        gwas_load_spec(),
        pd.DataFrame({"pos": [0], "p_value": [0.5]}),
        ["pos"],
        "'pos': 1 values <= 0",
    ),
    (
        "gwas-position-negative",
        gwas_load_spec(),
        pd.DataFrame({"pos": [-1], "p_value": [0.5]}),
        ["pos"],
        "'pos': 1 values <= 0",
    ),
    (
        "gwas-null-position",
        gwas_load_spec(),
        pd.DataFrame({"pos": [100, None], "p_value": [0.5, 0.5]}),
        ["pos"],
        "'pos' has 1 null values",
    ),
    (
        "gwas-null-pvalue",
        gwas_load_spec(),
        pd.DataFrame({"pos": [100, 200], "p_value": [0.05, None]}),
        ["p_value"],
        "'p_value' has 1 null values",
    ),
    (
        "gwas-pvalue-zero",
        gwas_load_spec(),
        pd.DataFrame({"pos": [100], "p_value": [0.0]}),
        ["p_value"],
        "'p_value': 1 values <= 0",
    ),
    (
        "gwas-pvalue-above-one",
        gwas_load_spec(),
        pd.DataFrame({"pos": [100], "p_value": [1.5]}),
        ["p_value"],
        "'p_value': 1 values > 1",
    ),
    (
        "gwas-non-numeric-pos",
        gwas_load_spec(),
        pd.DataFrame({"pos": ["chr1:1000", "chr1:2000"], "p_value": [0.01, 0.001]}),
        ["pos"],
        "must be numeric",
    ),
    (
        "gwas-non-numeric-pvalue",
        gwas_load_spec(),
        pd.DataFrame({"pos": [1000000, 1001000], "p_value": ["0.01", "significant"]}),
        ["p_value"],
        "must be numeric",
    ),
    (
        "gwas-both-non-numeric",
        gwas_load_spec(),
        pd.DataFrame({"pos": ["chr1:100"], "p_value": ["not_a_number"]}),
        ["pos", "p_value"],
        "must be numeric",
    ),
    (
        "eqtl-missing-all",
        EQTL_LOAD,
        pd.DataFrame({"extra": [1]}),
        ["pos", "p_value", "gene"],
        "Missing columns",
    ),
    (
        "eqtl-missing-gene-only",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100], "p_value": [0.05]}),
        ["gene"],
        "Missing columns",
    ),
    (
        "eqtl-non-numeric-pos",
        EQTL_LOAD,
        pd.DataFrame({"pos": ["chr1:100"], "p_value": [0.05], "gene": ["A"]}),
        ["pos"],
        "must be numeric",
    ),
    (
        "eqtl-non-numeric-pvalue",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100], "p_value": ["significant"], "gene": ["A"]}),
        ["p_value"],
        "must be numeric",
    ),
    (
        "eqtl-pvalue-zero",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100], "p_value": [0.0], "gene": ["A"]}),
        ["p_value"],
        "'p_value': 1 values <= 0",
    ),
    (
        "eqtl-non-positive-pos",
        EQTL_LOAD,
        pd.DataFrame({"pos": [-5, 100], "p_value": [0.05, 0.1], "gene": ["A", "B"]}),
        ["pos"],
        "'pos': 1 values <= 0",
    ),
    (
        "eqtl-null-pvalue",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100, 200], "p_value": [0.05, None], "gene": ["A", "B"]}),
        ["p_value"],
        "'p_value' has 1 null values",
    ),
    (
        "eqtl-null-pos",
        EQTL_LOAD,
        pd.DataFrame({"pos": [100, None], "p_value": [0.05, 0.1], "gene": ["A", "B"]}),
        ["pos"],
        "'pos' has 1 null values",
    ),
    (
        "finemapping-missing-both",
        FINEMAPPING_LOAD,
        pd.DataFrame({"other": [1]}),
        ["pos", "pip"],
        "Missing columns",
    ),
    (
        "finemapping-non-numeric-pip",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [100], "pip": ["high"]}),
        ["pip"],
        "must be numeric",
    ),
    (
        "finemapping-non-numeric-pos",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": ["abc"], "pip": [0.5]}),
        ["pos"],
        "must be numeric",
    ),
    (
        "finemapping-pip-above-one",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [100], "pip": [1.5]}),
        ["pip"],
        "'pip': 1 values > 1",
    ),
    (
        "finemapping-pip-negative",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [100], "pip": [-0.1]}),
        ["pip"],
        "'pip': 1 values < 0",
    ),
    (
        "finemapping-non-positive-pos",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [0], "pip": [0.5]}),
        ["pos"],
        "'pos': 1 values <= 0",
    ),
    (
        "finemapping-null-pip",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [100, 200], "pip": [0.5, None]}),
        ["pip"],
        "'pip' has 1 null values",
    ),
    (
        "finemapping-null-pos",
        FINEMAPPING_LOAD,
        pd.DataFrame({"pos": [100, None], "pip": [0.5, 0.6]}),
        ["pos"],
        "'pos' has 1 null values",
    ),
    (
        "genes-null-start",
        GENES_LOAD,
        _genes(
            chr=["1", "1"], start=[1000, None], end=[2000, 3000], gene_name=["A", "B"]
        ),
        ["start"],
        "'start' has 1 null values",
    ),
    (
        "genes-missing-all",
        GENES_LOAD,
        pd.DataFrame({"extra": [1]}),
        ["chr", "start", "end", "gene_name"],
        "Missing columns",
    ),
    (
        "genes-non-numeric-start",
        GENES_LOAD,
        _genes(start=["abc"]),
        ["start"],
        "must be numeric",
    ),
    (
        "genes-non-numeric-end",
        GENES_LOAD,
        _genes(end=["xyz"]),
        ["end"],
        "must be numeric",
    ),
    (
        "genes-negative-start",
        GENES_LOAD,
        _genes(start=[-1]),
        ["start"],
        "'start': 1 values < 0",
    ),
    (
        "genes-end-before-start",
        GENES_LOAD,
        _genes(start=[2000], end=[1000]),
        ["start", "end"],
        "rows have start > end",
    ),
]


@pytest.mark.parametrize(
    ("contract", "df", "named_columns", "fragment"),
    [
        pytest.param(f, df, cols, fragment, id=name)
        for name, f, df, cols, fragment in REJECTED
    ],
)
def test_rejected_input_raises_and_names_the_offending_columns(
    contract, df, named_columns, fragment
):
    """Rejection is a LoaderValidationError naming each bad column, and why."""
    with pytest.raises(LoaderValidationError) as exc_info:
        check(df, contract)

    message = str(exc_info.value)
    for column in named_columns:
        assert column in message, f"{column!r} missing from: {message}"
    if fragment is not None:
        assert fragment in message, f"{fragment!r} missing from: {message}"


def test_custom_gwas_column_names_are_honoured():
    """pos_col and p_col redirect validation to the named columns.

    A frame carrying the default names does not satisfy a spec asking for
    other ones, and the message names the column the spec asked for rather
    than the one the frame happens to have.
    """
    df = pd.DataFrame({"position": [100], "pvalue": [0.5]})
    check(df, gwas_load_spec(pos_col="position", p_col="pvalue"))

    with pytest.raises(LoaderValidationError) as exc_info:
        check(df, gwas_load_spec(pos_col="missing_pos", p_col="pvalue"))
    assert "missing_pos" in str(exc_info.value)

    defaults = pd.DataFrame({"pos": [100], "p_value": [0.5]})
    with pytest.raises(LoaderValidationError) as exc_info:
        check(defaults, gwas_load_spec(pos_col="position", p_col="p_value"))
    assert "Missing columns" in str(exc_info.value)
    assert "position" in str(exc_info.value)


def test_all_problems_are_reported_in_one_pass():
    """Validation accumulates rather than raising on the first fault."""
    df = pd.DataFrame({"pos": [0, 100], "p_value": [0.5, 2.0]})
    with pytest.raises(LoaderValidationError) as exc_info:
        check(df, gwas_load_spec())

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
            check(df, gwas_load_spec())
        except LoaderValidationError:
            return False
        return True

    @staticmethod
    def _kept_by_intake(value, *, family):
        from pylocuszoom._data import prepare_pvalue_data

        df = pd.DataFrame({"p": [value]})
        return len(prepare_pvalue_data(df, "p", family)) == 1

    @pytest.mark.parametrize(
        "value", [1e-300, 1e-8, 0.05, 0.5, 0.999, 1.0, 1.0000001, 1.5, 2.0, -0.1]
    )
    def test_strict_agrees_with_intake_when_zero_is_disallowed(self, value):
        """Away from zero, both ends of the library accept exactly the same values."""
        assert self._accepted_by_strict(value) == self._kept_by_intake(
            value, family="qq"
        )

    def test_zero_is_the_only_deliberate_divergence(self):
        """p == 0 is rejected at load and kept by the Manhattan convention."""
        assert self._accepted_by_strict(0.0) is False
        assert self._kept_by_intake(0.0, family="qq") is False
        assert self._kept_by_intake(0.0, family="genome-wide") is True

    def test_upper_bound_has_a_single_source(self):
        """The p-value spec rule and prepare_pvalue_data read the same constant."""
        from pylocuszoom._data import P_VALUE_MAX
        from pylocuszoom.exceptions import LoaderValidationError

        just_over = pd.DataFrame({"pos": [100], "p_value": [P_VALUE_MAX * 1.000001]})
        with pytest.raises(LoaderValidationError):
            check(just_over, gwas_load_spec())
        assert self._accepted_by_strict(P_VALUE_MAX) is True

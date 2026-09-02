"""Contract tests for the shared plot-data intake policy."""

import pandas as pd

from pylocuszoom._data import prepare_pvalue_data
from pylocuszoom.eqtl import prepare_eqtl_for_plotting
from pylocuszoom.manhattan import prepare_categorical_data, prepare_manhattan_data
from pylocuszoom.qq import prepare_qq_data

BAD_PVALUES = [1e-6, None, 0.0, -0.1, 1.5, 1e-320]


def _bad_pvalue_frame(p_col: str = "p") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["1"] * 6,
            "category": ["a"] * 6,
            "pos": [1, 2, 3, 4, 5, 6],
            p_col: BAD_PVALUES,
        }
    )


def test_eqtl_preparation_uses_the_shared_policy():
    raw = pd.DataFrame(
        {
            "pos": [1, 2, 3],
            "p_value": [1e-6, None, 2.0],
        }
    )
    result = prepare_eqtl_for_plotting(raw)
    expected = prepare_pvalue_data(raw, "p_value")
    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    ("entry_point", "allow_zero"),
    [
        (lambda df: prepare_manhattan_data(df, species="human"), True),
        (lambda df: prepare_categorical_data(df, "category"), True),
        (lambda df: transform_pvalues(df, "p"), True),
    ],
    ids=["manhattan", "categorical", "regional"],
)
def test_entry_points_keep_the_same_survivors(entry_point, allow_zero):
    """Every DataFrame-returning entry point drops the same bad rows."""
    raw = _bad_pvalue_frame()
    expected = prepare_pvalue_data(raw, "p", allow_zero=allow_zero)

    result = entry_point(raw)

    assert result["pos"].tolist() == expected["pos"].tolist()


def test_qq_excludes_exact_zero_pvalues():
    """QQ is the one entry point on the strict ``(0, 1]`` domain."""
    raw = _bad_pvalue_frame()
    expected = prepare_pvalue_data(raw, "p", allow_zero=False)

    result = prepare_qq_data(raw, p_col="p")

    assert result.attrs["n_variants"] == len(expected)


@pytest.mark.parametrize(
    ("entry_point", "message"),
    [
        (
            lambda df: prepare_manhattan_data(df, species="human"),
            "All rows have invalid p-values",
        ),
        (
            lambda df: prepare_categorical_data(df, "category"),
            "All rows have invalid p-values",
        ),
        (lambda df: prepare_qq_data(df, p_col="p"), "No valid p-values"),
    ],
    ids=["manhattan", "categorical", "qq"],
)
def test_entry_points_raise_when_nothing_survives(entry_point, message):
    """Manhattan-family and QQ raise rather than returning an empty frame."""
    raw = pd.DataFrame(
        {"chrom": ["1", "1"], "category": ["a", "a"], "pos": [1, 2], "p": [2.0, -1.0]}
    )

    with pytest.raises(ValueError, match=message):
        entry_point(raw)


def test_regional_returns_empty_when_nothing_survives():
    """The regional path returns an empty frame where the others raise."""
    raw = pd.DataFrame({"pos": [1, 2], "p": [2.0, -1.0]})

    assert transform_pvalues(raw, "p").empty

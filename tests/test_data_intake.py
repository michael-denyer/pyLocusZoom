"""Contract tests for the shared plot-data intake policy."""

import pandas as pd
import pytest

from pylocuszoom._data import prepare_pvalue_data
from pylocuszoom._plotter_utils import transform_pvalues
from pylocuszoom.eqtl import prepare_eqtl_for_plotting


def test_plot_families_share_pvalue_filtering_and_transform():
    raw = pd.DataFrame(
        {
            "pos": [1, 2, 3, 4, 5, 6],
            "p": [1e-6, None, 0.0, -0.1, 1.5, 1e-320],
        }
    )

    expected = prepare_pvalue_data(raw, "p")
    assert transform_pvalues(raw, "p")["pos"].tolist() == expected["pos"].tolist()
    assert transform_pvalues(raw, "p")["neglog10p"].tolist() == pytest.approx(
        expected["neglog10p"].tolist()
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

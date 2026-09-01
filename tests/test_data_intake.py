"""Contract tests for the shared plot-data intake policy."""

import pandas as pd

from pylocuszoom._data import prepare_pvalue_data
from pylocuszoom.eqtl import prepare_eqtl_for_plotting


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

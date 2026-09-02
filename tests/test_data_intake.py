"""Contract tests for the shared plot-data intake policy."""

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import settings as hyp_settings
from matplotlib.figure import Figure

from pylocuszoom._data import prepare_pvalue_data
from pylocuszoom.eqtl import prepare_eqtl_for_plotting
from pylocuszoom.manhattan import prepare_categorical_data, prepare_manhattan_data
from pylocuszoom.plotter import LocusZoomPlotter
from pylocuszoom.qq import prepare_qq_data
from tests.strategies import gwas_dataframes

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
    ],
    ids=["manhattan", "categorical"],
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


class TestPValueValidation:
    """Tests for p-value validation and NaN handling."""

    @pytest.fixture
    def speciesless_plotter(self):
        """Create plotter instance."""
        return LocusZoomPlotter(species=None)

    def test_plot_handles_nan_pvalues_with_warning(self, speciesless_plotter):
        """Plot should handle NaN p-values and log a warning."""
        import numpy as np

        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, np.nan, 0.01],  # One NaN p-value
            }
        )

        # Should not raise, but should warn (captured by logging)
        fig = speciesless_plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert isinstance(fig, Figure)

    def test_plot_stacked_handles_all_nan_pvalues(self, speciesless_plotter):
        """plot_stacked should handle region with all NaN p-values.

        Regression test: idxmin() on all-NaN series returns NaN,
        causing subsequent loc to fail.
        """
        import numpy as np

        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [np.nan, np.nan, np.nan],  # All NaN
            }
        )

        # Should not raise - should handle gracefully
        fig = speciesless_plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert isinstance(fig, Figure)

    def test_plot_handles_out_of_range_pvalues(self, speciesless_plotter):
        """Plot should handle p-values outside [0, 1] range."""
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [-0.1, 1.5, 0.05],  # Out of range values
            }
        )

        # Should not raise, but should warn
        fig = speciesless_plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert isinstance(fig, Figure)

    def test_prepare_pvalue_data_filters_invalid_range(self, speciesless_plotter):
        """prepare_pvalue_data filters out-of-range p-values (< 0 or > 1)."""
        import io

        from loguru import logger as loguru_logger

        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1700000, 1900000],
                "p_wald": [0.001, 0.5, -0.1, 1.5],
            }
        )

        # Capture log output
        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="WARNING",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        try:
            result = prepare_pvalue_data(df.copy(), "p_wald")
        finally:
            loguru_logger.remove(handler_id)

        # Only valid rows (0.001, 0.5) should remain
        assert len(result) == 2
        assert "neglog10p" in result.columns
        assert np.isfinite(result["neglog10p"]).all()
        assert (result["neglog10p"] > 0).all()

        # Warning logged about 2 invalid values
        log_output = log_capture.getvalue()
        assert "2 p-values outside [0, 1]" in log_output

    def test_prepare_pvalue_data_filters_nan(self, speciesless_plotter):
        """prepare_pvalue_data filters NaN p-values."""
        import io

        from loguru import logger as loguru_logger

        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [0.001, np.nan, 0.5],
            }
        )

        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="WARNING",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        try:
            result = prepare_pvalue_data(df.copy(), "p_wald")
        finally:
            loguru_logger.remove(handler_id)

        # NaN row filtered out
        assert len(result) == 2
        assert not result["p_wald"].isna().any()

        # Warning logged
        log_output = log_capture.getvalue()
        assert "1 NaN p-values" in log_output

    def test_prepare_pvalue_data_clips_very_small(self, speciesless_plotter):
        """prepare_pvalue_data clips very small p-values to 1e-300."""
        import io

        from loguru import logger as loguru_logger

        df = pd.DataFrame(
            {
                "ps": [1100000],
                "p_wald": [1e-310],  # Smaller than 1e-300
            }
        )

        log_capture = io.StringIO()
        handler_id = loguru_logger.add(
            log_capture,
            level="DEBUG",
            format="{message}",
            filter=lambda record: record["name"].startswith("pylocuszoom"),
        )

        try:
            result = prepare_pvalue_data(df.copy(), "p_wald")
        finally:
            loguru_logger.remove(handler_id)

        # Should be clipped to 1e-300, giving -log10(1e-300) = 300
        assert len(result) == 1
        assert result["neglog10p"].iloc[0] == pytest.approx(300.0)
        assert not np.isinf(result["neglog10p"]).any()

        # Debug log about clipping
        log_output = log_capture.getvalue()
        assert "Clipping" in log_output

    def test_prepare_pvalue_data_preserves_valid_data(self, speciesless_plotter):
        """prepare_pvalue_data passes through all-valid data unchanged."""
        df = pd.DataFrame(
            {
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [0.001, 0.5, 1e-8],
            }
        )

        result = prepare_pvalue_data(df.copy(), "p_wald")

        assert len(result) == 3
        assert result["neglog10p"].iloc[0] == pytest.approx(3.0)
        assert result["neglog10p"].iloc[1] == pytest.approx(-np.log10(0.5))
        assert result["neglog10p"].iloc[2] == pytest.approx(8.0)

    def test_plot_stacked_lead_detection_excludes_out_of_range(
        self, speciesless_plotter
    ):
        """Lead SNP auto-detection should exclude out-of-range p-values.

        Regression test: lead detection only filtered NaN but not out-of-range
        p-values, so a lead SNP with p=-0.1 could be selected and then removed
        by prepare_pvalue_data, causing missing lead highlighting.
        """
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "ps": [1100000, 1500000, 1900000],
                # rs1 has smallest absolute value but is invalid (negative)
                # rs3 should be selected as lead (smallest valid p-value)
                "p_wald": [-0.1, 0.5, 0.001],
            }
        )

        # Should not raise — lead detection should skip the invalid p-value
        fig = speciesless_plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert isinstance(fig, Figure)

    def test_plot_stacked_all_invalid_pvalues(self, speciesless_plotter):
        """plot_stacked should handle region with all out-of-range p-values."""
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "chr": [1, 1],
                "ps": [1100000, 1500000],
                "p_wald": [-0.1, 1.5],  # All invalid
            }
        )

        fig = speciesless_plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert isinstance(fig, Figure)


class TestPvalueTransformation:
    """Tests for p-value transformation helper."""

    def test_prepare_pvalue_data_adds_neglog10p_column(self):
        """Helper creates neglog10p column from p-values."""
        df = pd.DataFrame({"pval": [0.01, 0.001, 1e-8]})

        result = prepare_pvalue_data(df.copy(), "pval")

        assert "neglog10p" in result.columns
        assert result["neglog10p"].iloc[0] == pytest.approx(2.0)  # -log10(0.01)
        assert result["neglog10p"].iloc[1] == pytest.approx(3.0)  # -log10(0.001)
        assert result["neglog10p"].iloc[2] == pytest.approx(8.0)  # -log10(1e-8)

    def test_prepare_pvalue_data_clips_extreme_values(self):
        """Extremely small p-values are clipped to avoid -inf."""
        df = pd.DataFrame({"pval": [1e-350, 0.0]})  # Would be -inf without clipping

        result = prepare_pvalue_data(df.copy(), "pval")

        # Should be clipped to 1e-300, giving ~300
        assert result["neglog10p"].iloc[0] == pytest.approx(300.0)
        assert result["neglog10p"].iloc[1] == pytest.approx(300.0)
        assert not np.isinf(result["neglog10p"]).any()


class TestPlotStackedProperties:
    """Property-based tests for stacked plots."""

    @hyp_settings(max_examples=10, deadline=None)
    @given(gwas_dataframes(min_snps=5, max_snps=30))
    def test_plot_stacked_with_duplicated_data(self, df):
        """plot_stacked() should handle multiple identical DataFrames."""
        plotter = LocusZoomPlotter(species="canine")
        chrom = df["chr"].iloc[0]
        start = int(df["ps"].min())
        end = int(df["ps"].max())

        fig = plotter.plot_stacked(
            [df, df.copy()],
            chrom=chrom,
            start=start,
            end=end,
            show_recombination=False,
        )

        assert fig is not None

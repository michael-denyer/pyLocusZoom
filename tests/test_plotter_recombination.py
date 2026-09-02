"""Tests for the recombination overlay in regional plots."""

from unittest.mock import patch

import pandas as pd
import pytest

from pylocuszoom.plotter import LocusZoomPlotter
from pylocuszoom.recombination import RecombResult, RecombStatus


class TestLocusZoomPlotterRecombination:
    """Tests for recombination data handling."""

    def test_caches_recombination_data(self):
        """Should cache recombination data for repeated calls."""
        plotter = LocusZoomPlotter(species=None)  # No auto-download

        recomb_df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "rate": [0.5, 1.0, 0.5],
            }
        )

        # First call - no cache
        assert plotter._recomb_cache == {}

        # Manually add to cache (key includes genome_build)
        plotter._recomb_cache[(1, 1000000, 2000000, plotter.genome_build)] = (
            RecombResult(RecombStatus.OK, frame=recomb_df)
        )

        # Should return cached data
        result = plotter._get_recomb_for_region(1, 1000000, 2000000)
        assert result.status is RecombStatus.OK
        assert len(result.frame) == 3

    def test_recombination_overlay_does_not_distort_primary_ylim(self):
        """Primary y-axis limits should be unchanged when recombination is enabled.

        Regression test: recombination overlay was being plotted on the primary axis
        instead of a twin axis, causing GWAS y-limits to be rescaled by recomb rates.
        """
        plotter = LocusZoomPlotter(species=None)

        gwas_df = pd.DataFrame(
            {
                "rs": [f"rs{i}" for i in range(10)],
                "chr": [1] * 10,
                "ps": list(range(1000000, 2000000, 100000)),
                "p_wald": [1e-8, 1e-6, 1e-5, 1e-4, 0.01, 0.05, 0.1, 0.5, 0.8, 0.99],
            }
        )

        recomb_df = pd.DataFrame(
            {
                "pos": [1000000, 1500000, 2000000],
                "rate": [50.0, 100.0, 75.0],  # High rates that would distort y-axis
            }
        )

        # Plot without recombination
        fig_no_recomb = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        ax_no_recomb = fig_no_recomb.axes[0]
        ylim_no_recomb = ax_no_recomb.get_ylim()

        # Plot with recombination
        fig_with_recomb = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            recomb_df=recomb_df,
        )
        ax_with_recomb = fig_with_recomb.axes[0]
        ylim_with_recomb = ax_with_recomb.get_ylim()

        # Primary y-axis limits should be the same
        assert ylim_no_recomb == ylim_with_recomb, (
            f"Recombination overlay distorted primary y-axis: "
            f"without={ylim_no_recomb}, with={ylim_with_recomb}"
        )


class TestRecombinationDownloadErrors:
    """Tests for recombination map error handling.

    These tests verify that when recombination maps are unavailable,
    the plotter gracefully handles None return values and allows
    plotting to continue without recombination overlay.

    Note: Detailed error handling (network, I/O, OS errors) is tested
    in test_recombination.py at the recomb_for_region level.
    """

    @pytest.fixture
    def debug_canine_plotter(self):
        """Create a plotter instance for testing download errors."""
        return LocusZoomPlotter(species="canine", log_level="DEBUG")

    @staticmethod
    def _download_failed():
        return patch(
            "pylocuszoom.plotter.recomb_for_region",
            return_value=RecombResult(
                RecombStatus.DOWNLOAD_FAILED, detail="could not download maps"
            ),
        )

    def test_plotting_continues_without_recomb_maps(
        self, debug_canine_plotter, tiny_regional_gwas_df
    ):
        """Plotting should succeed even when recombination maps are unavailable."""
        with self._download_failed():
            with pytest.warns(UserWarning, match="could not download maps"):
                fig = debug_canine_plotter.plot(
                    tiny_regional_gwas_df,
                    chrom=1,
                    start=1000000,
                    end=2000000,
                    show_recombination=True,
                )
            assert fig is not None

    def test_a_skipped_overlay_warns_once(
        self, debug_canine_plotter, tiny_regional_gwas_df
    ):
        """Three layers used to decide this; only one of them speaks now."""
        with self._download_failed(), pytest.warns(UserWarning) as caught:
            debug_canine_plotter.plot(
                tiny_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=True,
            )

        skipped = [
            w for w in caught if "Recombination overlay skipped" in str(w.message)
        ]
        assert len(skipped) == 1

    def test_the_warning_points_at_the_caller_not_the_library(
        self, debug_canine_plotter, tiny_regional_gwas_df
    ):
        """A file:line inside pylocuszoom tells the user nothing actionable."""
        with self._download_failed(), pytest.warns(UserWarning) as caught:
            debug_canine_plotter.plot(
                tiny_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=True,
            )

        skipped = next(
            w for w in caught if "Recombination overlay skipped" in str(w.message)
        )
        assert skipped.filename == __file__


class TestRecombinationOptionalDependency:
    """The overlay is skipped when pyliftover is missing and propagates any
    other ImportError, decided by exception type rather than message text."""

    @pytest.fixture
    def plotter(self, tmp_path):
        return LocusZoomPlotter(species="canine", log_level=None)

    @pytest.fixture(autouse=True)
    def _maps_are_present(self, tmp_path):
        with patch(
            "pylocuszoom.recombination.ensure_recomb_maps", return_value=tmp_path
        ):
            yield

    def test_missing_optional_dependency_is_reported_as_a_status(self, plotter):
        from pylocuszoom.exceptions import OptionalDependencyMissing

        with patch(
            "pylocuszoom.recombination.get_recombination_rate_for_region",
            side_effect=OptionalDependencyMissing("no liftover here"),
        ):
            result = plotter._get_recomb_for_region(1, 1_000_000, 2_000_000)

        assert result.status is RecombStatus.LIFTOVER_UNAVAILABLE
        assert "no liftover here" in result.detail

    def test_other_import_error_propagates(self, plotter):
        with patch(
            "pylocuszoom.recombination.get_recombination_rate_for_region",
            side_effect=ImportError("pyliftover mentioned but unrelated"),
        ):
            with pytest.raises(ImportError):
                plotter._get_recomb_for_region(1, 1_000_000, 2_000_000)

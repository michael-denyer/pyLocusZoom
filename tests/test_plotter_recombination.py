"""Tests for the recombination overlay in regional plots."""

from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pylocuszoom.plotter import LocusZoomPlotter


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
        plotter._recomb_cache[(1, 1000000, 2000000, plotter.genome_build)] = recomb_df

        # Should return cached data
        result = plotter._get_recomb_for_region(1, 1000000, 2000000)
        assert result is not None
        assert len(result) == 3

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

        plt.close(fig_no_recomb)
        plt.close(fig_with_recomb)


class TestRecombinationDownloadErrors:
    """Tests for recombination map error handling.

    These tests verify that when recombination maps are unavailable,
    the plotter gracefully handles None return values and allows
    plotting to continue without recombination overlay.

    Note: Detailed error handling (network, I/O, OS errors) is tested
    in test_recombination.py at the ensure_recomb_maps level.
    """

    @pytest.fixture
    def debug_canine_plotter(self):
        """Create a plotter instance for testing download errors."""
        return LocusZoomPlotter(species="canine", log_level="DEBUG")

    def test_ensure_recomb_maps_returns_none_propagates(self, debug_canine_plotter):
        """When ensure_recomb_maps returns None, the plotter returns None."""
        with patch(
            "pylocuszoom.plotter.ensure_recomb_maps", return_value=None
        ) as mock_ensure:
            result = debug_canine_plotter._ensure_recomb_maps()
            assert result is None
            mock_ensure.assert_called_once_with(species="canine", data_dir=None)

    def test_plotting_continues_without_recomb_maps(
        self, debug_canine_plotter, tiny_regional_gwas_df
    ):
        """Plotting should succeed even when recombination maps are unavailable."""
        with patch("pylocuszoom.plotter.ensure_recomb_maps", return_value=None):
            # Should not raise, just skip recombination overlay
            fig = debug_canine_plotter.plot(
                tiny_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=True,
            )
            assert fig is not None
            plt.close(fig)

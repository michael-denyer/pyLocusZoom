"""Tests for the recombination overlay in regional plots."""

from unittest.mock import patch

import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, PanelInputs
from pylocuszoom.plotter import LocusZoomPlotter
from pylocuszoom.recombination import RecombResult, RecombStatus


class TestLocusZoomPlotterRecombination:
    """Tests for recombination data handling."""

    def test_caches_recombination_data(self, canine_plotter, tmp_path):
        """A region plotted twice reads its maps once."""
        region = dict(
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(snp_labels=False),
        )
        gwas_df = pd.DataFrame(
            {"chr": 1, "pos": [1100000, 1900000], "p_value": [1e-8, 1e-3]}
        )

        def recomb_rates(fig):
            (line,) = fig.axes[1].get_lines()
            return list(line.get_ydata())

        first = recomb_rates(canine_plotter.plot(gwas_df, **region))
        (tmp_path / "recombination_maps" / "chr1_recomb.tsv").write_text(
            "chr\tpos\trate\tcM\n1\t1500000\t9.0\t1.0\n"
        )
        again = recomb_rates(canine_plotter.plot(gwas_df, **region))
        fresh = LocusZoomPlotter(
            species="canine", recomb_data_dir=canine_plotter.recomb_data_dir
        )

        assert first == again == [1.0]
        assert recomb_rates(fresh.plot(gwas_df, **region)) == [9.0]

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
                "pos": list(range(1000000, 2000000, 100000)),
                "p_value": [1e-8, 1e-6, 1e-5, 1e-4, 0.01, 0.05, 0.1, 0.5, 0.8, 0.99],
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
            display=DisplayConfig(show_recombination=False),
        )
        ax_no_recomb = fig_no_recomb.axes[0]
        ylim_no_recomb = ax_no_recomb.get_ylim()

        # Plot with recombination
        fig_with_recomb = plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            panels=PanelInputs(recomb_df=recomb_df),
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
                    display=DisplayConfig(show_recombination=True),
                )
        assert [ax.get_ylabel() for ax in fig.axes] == [r"$-\log_{10}$ P"]

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
                display=DisplayConfig(show_recombination=True),
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
                display=DisplayConfig(show_recombination=True),
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

    def test_missing_optional_dependency_skips_the_overlay_with_a_warning(
        self, plotter, tiny_regional_gwas_df
    ):
        from pylocuszoom.exceptions import OptionalDependencyMissing

        with (
            patch(
                "pylocuszoom.recombination.get_recombination_rate_for_region",
                side_effect=OptionalDependencyMissing("no liftover here"),
            ),
            pytest.warns(UserWarning, match="no liftover here"),
        ):
            fig = plotter.plot(
                tiny_regional_gwas_df, chrom=1, start=1_000_000, end=2_000_000
            )

        assert [ax.get_ylabel() for ax in fig.axes] == [r"$-\log_{10}$ P"]

    def test_other_import_error_propagates(self, plotter, tiny_regional_gwas_df):
        with patch(
            "pylocuszoom.recombination.get_recombination_rate_for_region",
            side_effect=ImportError("pyliftover mentioned but unrelated"),
        ):
            with pytest.raises(ImportError, match="mentioned but unrelated"):
                plotter.plot(
                    tiny_regional_gwas_df, chrom=1, start=1_000_000, end=2_000_000
                )


def test_a_failed_chain_download_warns_and_still_plots(
    cache_home, monkeypatch, tiny_regional_gwas_df
):
    """The chain is part of the overlay; losing it must not lose the figure."""
    from pylocuszoom.exceptions import DataDownloadError
    from tests.conftest import write_canine_map_set

    write_canine_map_set(
        cache_home / "recombination_maps",
        "chr\tpos\trate\tcM\n1\t1500000\t1.0\t0.1\n",
    )

    def refuse(*args, **kwargs):
        raise DataDownloadError("simulated chain 404")

    monkeypatch.setattr("pylocuszoom._liftover.download_file", refuse)
    plotter = LocusZoomPlotter(species="canine", genome_build="canfam4", log_level=None)

    with pytest.warns(UserWarning, match="simulated chain 404"):
        fig = plotter.plot(tiny_regional_gwas_df, chrom=1, start=1000000, end=2000000)

    assert fig.get_axes()

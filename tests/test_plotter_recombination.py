"""Tests for the recombination overlay in regional plots."""

from unittest.mock import patch

import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, PanelInputs
from pylocuszoom.exceptions import DataDownloadError
from pylocuszoom.plotter import LocusZoomPlotter
from tests.conftest import write_canine_map_set


def _skip_warnings(caught):
    return [w for w in caught if "Recombination overlay skipped" in str(w.message)]


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
    """A map set that cannot be had skips the overlay, never the figure."""

    @pytest.fixture
    def canine(self, cache_home):
        return LocusZoomPlotter(species="canine", log_level=None)

    @staticmethod
    def _download_failed():
        return patch(
            "pylocuszoom.recombination.download_file",
            side_effect=DataDownloadError("could not download maps"),
        )

    def test_plotting_continues_without_recomb_maps(
        self, canine, tiny_regional_gwas_df
    ):
        with self._download_failed():
            with pytest.warns(UserWarning, match="could not download maps"):
                fig = canine.plot(
                    tiny_regional_gwas_df, chrom=1, start=1000000, end=2000000
                )
        assert [ax.get_ylabel() for ax in fig.axes] == [r"$-\log_{10}$ P"]

    def test_a_skipped_overlay_warns_once_and_points_at_the_caller(
        self, canine, tiny_regional_gwas_df
    ):
        """A file:line inside pylocuszoom tells the user nothing actionable."""
        with self._download_failed(), pytest.warns(UserWarning) as caught:
            canine.plot(tiny_regional_gwas_df, chrom=1, start=1000000, end=2000000)

        assert len(caught) == 1
        assert len(_skip_warnings(caught)) == 1
        assert caught[0].filename == __file__

    def test_plot_stacked_warns_at_its_caller_too(self, canine, tiny_regional_gwas_df):
        with self._download_failed(), pytest.warns(UserWarning) as caught:
            canine.plot_stacked(
                [tiny_regional_gwas_df], chrom=1, start=1000000, end=2000000
            )

        assert [w.filename for w in caught] == [__file__]

    def test_a_failed_download_is_retried_on_the_next_plot(
        self, canine, cache_home, tiny_regional_gwas_df
    ):
        """Failures are not memoised; a region that failed once can recover."""
        with self._download_failed(), pytest.warns(UserWarning):
            canine.plot(tiny_regional_gwas_df, chrom=1, start=1000000, end=2000000)
        write_canine_map_set(
            cache_home / "recombination_maps",
            "chr\tpos\trate\tcM\n1\t1500000\t3.0\t0.1\n",
        )

        fig = canine.plot(tiny_regional_gwas_df, chrom=1, start=1000000, end=2000000)

        assert fig.axes[1].get_ylabel() == "Recombination rate (cM/Mb)"

    def test_an_unwritable_cache_skips_the_overlay_with_one_warning(
        self, tmp_path, monkeypatch, tiny_regional_gwas_df
    ):
        blocker = tmp_path / "not_a_directory"
        blocker.write_text("")
        monkeypatch.setenv("XDG_CACHE_HOME", str(blocker / "cache"))
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        with pytest.warns(UserWarning) as caught:
            fig = plotter.plot(
                tiny_regional_gwas_df, chrom=1, start=1_000_000, end=2_000_000
            )

        assert len(caught) == 1
        assert "Could not write recombination maps" in str(caught[0].message)
        assert fig.get_axes()


class TestRecombinationOptionalDependency:
    """The overlay is skipped when pyliftover is missing and propagates any
    other ImportError, decided by exception type rather than message text."""

    @pytest.fixture
    def plotter(self, cache_home, monkeypatch):
        write_canine_map_set(
            cache_home / "recombination_maps",
            "chr\tpos\trate\tcM\n1\t1500000\t1.0\t0.1\n",
        )
        return LocusZoomPlotter(
            species="canine", genome_build="canfam4", log_level=None
        )

    def test_missing_optional_dependency_skips_the_overlay_with_a_warning(
        self, plotter, monkeypatch, tiny_regional_gwas_df
    ):
        import gzip
        import sys

        def download(url, dest, desc=None):
            dest.write_bytes(gzip.compress(b"chain"))

        monkeypatch.setattr("pylocuszoom._liftover.download_file", download)
        monkeypatch.setitem(sys.modules, "pyliftover", None)

        with pytest.warns(UserWarning, match="pip install pyliftover") as caught:
            fig = plotter.plot(
                tiny_regional_gwas_df, chrom=1, start=1_000_000, end=2_000_000
            )

        assert len(caught) == 1
        assert [ax.get_ylabel() for ax in fig.axes] == [r"$-\log_{10}$ P"]

    def test_other_import_error_propagates(
        self, plotter, monkeypatch, tiny_regional_gwas_df
    ):
        def broken(*args, **kwargs):
            raise ImportError("pyliftover mentioned but unrelated")

        monkeypatch.setattr("pylocuszoom._liftover.download_file", broken)

        with pytest.raises(ImportError, match="mentioned but unrelated"):
            plotter.plot(tiny_regional_gwas_df, chrom=1, start=1_000_000, end=2_000_000)


def test_a_failed_chain_download_warns_once_and_still_plots(
    cache_home, monkeypatch, tiny_regional_gwas_df
):
    """The chain is part of the overlay; losing it must not lose the figure."""
    write_canine_map_set(
        cache_home / "recombination_maps",
        "chr\tpos\trate\tcM\n1\t1500000\t1.0\t0.1\n",
    )

    def refuse(*args, **kwargs):
        raise DataDownloadError("simulated chain 404")

    monkeypatch.setattr("pylocuszoom._liftover.download_file", refuse)
    plotter = LocusZoomPlotter(species="canine", genome_build="canfam4", log_level=None)

    with pytest.warns(UserWarning, match="simulated chain 404") as caught:
        fig = plotter.plot(tiny_regional_gwas_df, chrom=1, start=1000000, end=2000000)

    assert len(caught) == 1
    assert fig.get_axes()

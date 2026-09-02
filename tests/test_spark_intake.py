"""Every public plot method collects its frames through ``utils.to_pandas``.

``docs/ARCHITECTURE.md`` describes the intake step as normalising through
``to_pandas``. These tests are what make that true: each plotter is handed an
object that is not a pandas DataFrame and only answers ``toPandas()``, and the
figure still comes back.
"""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom import (
    ColocPlotter,
    LocusZoomPlotter,
    ManhattanPlotter,
    MiamiPlotter,
    StatsPlotter,
)


class FakeSparkFrame:
    """A frame-like object with the one method the Spark contract needs."""

    def __init__(self, frame: pd.DataFrame):
        self._frame = frame
        self.collected = 0

    def toPandas(self) -> pd.DataFrame:  # noqa: N802 - the Spark spelling
        self.collected += 1
        return self._frame


@pytest.fixture
def spark_regional(tiny_regional_gwas_df):
    return FakeSparkFrame(tiny_regional_gwas_df)


@pytest.fixture
def spark_manhattan(manhattan_gwas_df):
    return FakeSparkFrame(manhattan_gwas_df)


def test_regional_plot_collects_the_frame(spark_regional):
    """plot() reaches the association panel from a Spark-like frame."""
    plotter = LocusZoomPlotter(species=None, log_level=None)

    fig = plotter.plot(spark_regional, chrom=1, start=1_000_000, end=2_000_000)

    assert fig is not None
    assert spark_regional.collected == 1


def test_stacked_plot_collects_every_frame(tiny_regional_gwas_df):
    """plot_stacked() collects each panel's frame, not only the first."""
    frames = [FakeSparkFrame(tiny_regional_gwas_df) for _ in range(2)]
    plotter = LocusZoomPlotter(species=None, log_level=None)

    plotter.plot_stacked(frames, chrom=1, start=1_000_000, end=2_000_000)

    assert [f.collected for f in frames] == [1, 1]


def test_manhattan_collects_the_frame(spark_manhattan):
    """plot_manhattan() lays out a Spark-like frame."""
    assert ManhattanPlotter(species="human").plot_manhattan(spark_manhattan)
    assert spark_manhattan.collected == 1


def test_qq_collects_the_frame(spark_manhattan):
    """plot_qq() reads p-values from a Spark-like frame."""
    assert ManhattanPlotter(species="human").plot_qq(spark_manhattan)


def test_manhattan_stacked_collects_every_frame(manhattan_gwas_df):
    """The stacked genome-wide entry points collect each frame."""
    frames = [FakeSparkFrame(manhattan_gwas_df) for _ in range(2)]

    ManhattanPlotter(species="human").plot_manhattan_stacked(frames)

    assert [f.collected for f in frames] == [1, 1]


def test_miami_collects_both_halves(manhattan_gwas_df):
    """plot_miami() collects the top and bottom frames."""
    top, bottom = (FakeSparkFrame(manhattan_gwas_df) for _ in range(2))

    MiamiPlotter(species="human").plot_miami(top, bottom)

    assert (top.collected, bottom.collected) == (1, 1)


def test_phewas_collects_the_frame(phewas_with_effects_df):
    """plot_phewas() reads a Spark-like frame."""
    frame = FakeSparkFrame(phewas_with_effects_df)

    StatsPlotter().plot_phewas(frame, variant_id="rs1", p_col="p")

    assert frame.collected == 1


def test_forest_collects_the_frame(sample_forest_df):
    """plot_forest() reads a Spark-like frame."""
    frame = FakeSparkFrame(sample_forest_df)

    StatsPlotter().plot_forest(frame, variant_id="rs1")

    assert frame.collected == 1


def test_coloc_collects_both_frames():
    """plot_coloc() collects the GWAS and eQTL frames."""
    rng = np.random.default_rng(0)
    positions = np.arange(1_000_000, 1_000_050)
    gwas = FakeSparkFrame(
        pd.DataFrame({"pos": positions, "p_gwas": rng.uniform(1e-9, 1, 50)})
    )
    eqtl = FakeSparkFrame(
        pd.DataFrame({"pos": positions, "p_eqtl": rng.uniform(1e-9, 1, 50)})
    )

    ColocPlotter().plot_coloc(gwas, eqtl, rs_col=None)

    assert (gwas.collected, eqtl.collected) == (1, 1)

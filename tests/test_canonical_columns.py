"""The canonical column vocabulary and its one-release deprecation path."""

import warnings

import pandas as pd
import pytest

from pylocuszoom import (
    ColumnConfig,
    GenomeWideConfig,
    LocusZoomPlotter,
    ManhattanPlotter,
    load_gemma,
)
from pylocuszoom.exceptions import ValidationError
from pylocuszoom.schemas import (
    DEPRECATED_ALIAS_REMOVED_IN,
    DEPRECATED_COLUMN_ALIASES,
    Canonical,
)

GEMMA_FILE = """chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tp_wald
1\trs1\t1100000\t0\tA\tG\t0.3\t0.5\t0.2\t1e-8
1\trs2\t1500000\t0\tA\tG\t0.3\t0.4\t0.2\t1e-5
1\trs3\t1900000\t0\tA\tG\t0.3\t0.3\t0.2\t1e-3
"""


@pytest.fixture
def gemma_file(tmp_path):
    """A three-variant GEMMA file in the format's own column names."""
    path = tmp_path / "output.assoc.txt"
    path.write_text(GEMMA_FILE)
    return path


@pytest.fixture
def legacy_gwas_df():
    """A region frame in the pre-4.0 loader output names."""
    return pd.DataFrame(
        {
            "chr": [1, 1, 1],
            "rs": ["rs1", "rs2", "rs3"],
            "ps": [1100000, 1500000, 1900000],
            "p_wald": [1e-8, 1e-5, 1e-3],
        }
    )


class TestCanonicalVocabulary:
    """One spelling per concept, across the loaders and the config models."""

    def test_loaders_emit_the_canonical_names(self, gemma_file):
        """A GWAS file loads into chr/pos/p_value/rs whatever the format calls them."""
        df = load_gemma(gemma_file)

        assert {Canonical.CHROM, Canonical.POS, Canonical.P, Canonical.RS} <= set(
            df.columns
        )
        assert "ps" not in df.columns
        assert "p_wald" not in df.columns

    def test_loading_emits_no_warning(self, gemma_file):
        """The canonical path is the quiet path."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            load_gemma(gemma_file)

    def test_the_config_models_default_to_it(self):
        """Both column models name the same columns the loaders emit."""
        assert ColumnConfig().pos_col == Canonical.POS
        assert ColumnConfig().p_col == Canonical.P
        assert ColumnConfig().rs_col == Canonical.RS
        assert GenomeWideConfig().chrom_col == Canonical.CHROM
        assert GenomeWideConfig().pos_col == Canonical.POS
        assert GenomeWideConfig().p_col == Canonical.P

    def test_a_loaded_frame_plots_without_renaming(self, gemma_file):
        """Draw every loaded variant at the position and p-value the loader emitted."""
        df = load_gemma(gemma_file)

        fig = ManhattanPlotter(species="human").plot_manhattan(df)

        ax = fig.get_axes()[0]
        drawn = {
            (float(x), round(float(y), 6))
            for collection in ax.collections
            for x, y in collection.get_offsets()
        }
        assert drawn == {(1100000.0, 8.0), (1500000.0, 5.0), (1900000.0, 3.0)}


class TestDeprecatedLoaderKnobs:
    """Naming a loader's output columns works for one more release."""

    def test_naming_a_column_still_renames_it(self, gemma_file):
        """The knob keeps its old behaviour while it is deprecated."""
        with pytest.warns(DeprecationWarning):
            df = load_gemma(gemma_file, pos_col="position", p_col="pvalue")

        assert "position" in df.columns
        assert "pvalue" in df.columns

    def test_the_warning_names_the_removal_release(self, gemma_file):
        """A user needs to know how long they have."""
        with pytest.warns(DeprecationWarning, match=DEPRECATED_ALIAS_REMOVED_IN):
            load_gemma(gemma_file, pos_col="position")

    def test_the_warning_names_the_knob(self, gemma_file):
        """Only the knobs actually passed are named."""
        with pytest.warns(DeprecationWarning, match=r"^p_col on the GWAS loaders"):
            load_gemma(gemma_file, p_col="pvalue")


class TestDeprecatedFrameColumns:
    """A frame from a 3.x loader still plots, once, with a warning."""

    def test_the_alias_table_covers_position_and_pvalue(self):
        """Chromosome and SNP id never changed spelling, so they have no alias."""
        assert DEPRECATED_COLUMN_ALIASES == {
            Canonical.POS: "ps",
            Canonical.P: "p_wald",
        }

    def test_regional_plot_accepts_the_old_names(self, legacy_gwas_df):
        """plot() falls back to ps/p_wald and says so."""
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with pytest.warns(DeprecationWarning, match="p_wald"):
            fig = plotter.plot(legacy_gwas_df, chrom=1, start=1_000_000, end=2_000_000)

        assert fig is not None

    def test_the_frame_warning_names_the_canonical_column(self, legacy_gwas_df):
        """The message tells the user what to rename to."""
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with pytest.warns(DeprecationWarning, match=r"pre-4\.0 name for 'pos'"):
            plotter.plot(legacy_gwas_df, chrom=1, start=1_000_000, end=2_000_000)

    def test_manhattan_accepts_the_old_names(self, legacy_gwas_df):
        """The genome-wide boundary takes the same fallback."""
        with pytest.warns(DeprecationWarning, match="ps"):
            fig = ManhattanPlotter(species="human").plot_manhattan(legacy_gwas_df)

        assert fig is not None

    def test_a_canonical_frame_does_not_warn(self, gemma_file):
        """The fallback fires on the old names only."""
        df = load_gemma(gemma_file)
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            plotter.plot(df, chrom=1, start=1_000_000, end=2_000_000)

    def test_a_canonical_column_wins_over_its_alias(self, legacy_gwas_df):
        """A frame carrying both is read as canonical, silently."""
        both = legacy_gwas_df.assign(pos=legacy_gwas_df["ps"], p_value=1e-9)
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            plotter.plot(both, chrom=1, start=1_000_000, end=2_000_000)

    def test_a_caller_named_column_gets_no_fallback(self, legacy_gwas_df):
        """Only the canonical names have aliases, so an explicit name is honoured."""
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with pytest.raises(ValidationError, match="position"):
            plotter.plot(
                legacy_gwas_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                columns=ColumnConfig(pos_col="position", p_col="p_wald"),
            )

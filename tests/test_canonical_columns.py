"""The canonical column vocabulary, and the 3.x spellings it no longer aliases."""

import warnings

import pandas as pd
import pytest

from pylocuszoom import (
    ColumnConfig,
    DisplayConfig,
    GenomeWideConfig,
    LocusZoomPlotter,
    ManhattanPlotter,
    load_gemma,
    load_gwas,
)
from pylocuszoom.exceptions import ValidationError
from pylocuszoom.schemas import Canonical

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


class TestRemovedLoaderKnobs:
    """The GWAS loaders' output-column overrides were removed in 5.0."""

    @pytest.mark.parametrize("knob", ["pos_col", "p_col", "rs_col"])
    def test_naming_an_output_column_is_an_error(self, gemma_file, knob):
        """A loader emits the canonical names; rename after loading instead."""
        with pytest.raises(TypeError, match=knob):
            load_gemma(gemma_file, **{knob: "renamed"})

    def test_load_gwas_takes_no_extra_arguments(self, gemma_file):
        """load_gwas(**kwargs) promised pass-through that no loader accepted."""
        with pytest.raises(TypeError, match="sep"):
            load_gwas(gemma_file, sep=",")


class TestPre4ColumnNamesAreNotAliases:
    """A frame in the 3.x ps/p_wald names is an ordinary frame missing columns."""

    def test_regional_plot_names_the_missing_canonical_column(self, legacy_gwas_df):
        plotter = LocusZoomPlotter(species=None)

        with pytest.raises(ValidationError, match="'pos', 'p_value'"):
            plotter.plot(legacy_gwas_df, chrom=1, start=1_000_000, end=2_000_000)

    def test_manhattan_names_the_missing_canonical_column(self, legacy_gwas_df):
        with pytest.raises(ValidationError, match="'pos', 'p_value'"):
            ManhattanPlotter(species="human").plot_manhattan(legacy_gwas_df)

    def test_qq_names_the_missing_canonical_column(self, legacy_gwas_df):
        with pytest.raises(ValidationError, match="p_value"):
            ManhattanPlotter(species="human").plot_qq(legacy_gwas_df)

    def test_naming_the_old_columns_plots_them(self, legacy_gwas_df):
        """The migration is one ColumnConfig, or one rename."""
        fig = LocusZoomPlotter(species=None).plot(
            legacy_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
            display=DisplayConfig(show_recombination=False),
        )

        heights = {
            round(float(y), 6)
            for collection in fig.axes[0].collections
            for _, y in collection.get_offsets()
        }
        assert heights == {8.0, 5.0, 3.0}

    def test_a_caller_named_column_gets_no_fallback(self, legacy_gwas_df):
        """An explicit name is honoured, not swapped for a similar column."""
        plotter = LocusZoomPlotter(species=None)

        with pytest.raises(ValidationError, match="position"):
            plotter.plot(
                legacy_gwas_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                columns=ColumnConfig(pos_col="position", p_col="p_wald"),
            )

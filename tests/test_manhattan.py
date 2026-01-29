"""Tests for Manhattan plot functionality."""

import numpy as np
import pandas as pd
import pytest


class TestChromosomeOrdering:
    """Tests for chromosome order handling."""

    def test_canine_chromosome_order_includes_autosomes(self):
        """Canine order should have 38 autosomes."""
        from pylocuszoom.manhattan import get_chromosome_order

        order = get_chromosome_order(species="canine")
        autosomes = [str(i) for i in range(1, 39)]
        for chrom in autosomes:
            assert chrom in order, f"Missing autosome {chrom}"

    def test_canine_chromosome_order_includes_sex_and_mt(self):
        """Canine order should include X, Y, MT."""
        from pylocuszoom.manhattan import get_chromosome_order

        order = get_chromosome_order(species="canine")
        assert "X" in order
        assert "Y" in order
        assert "MT" in order

    def test_feline_chromosome_order_uses_letter_format(self):
        """Feline chromosomes use A1-E3 format."""
        from pylocuszoom.manhattan import get_chromosome_order

        order = get_chromosome_order(species="feline")
        assert "A1" in order
        assert "B4" in order
        assert "E3" in order
        assert "X" in order
        assert "Y" in order
        assert "MT" in order

    def test_human_chromosome_order_has_22_autosomes(self):
        """Human order should have 22 autosomes."""
        from pylocuszoom.manhattan import get_chromosome_order

        order = get_chromosome_order(species="human")
        autosomes = [str(i) for i in range(1, 23)]
        for chrom in autosomes:
            assert chrom in order

    def test_custom_order_overrides_species(self):
        """Custom order should be used when provided."""
        from pylocuszoom.manhattan import get_chromosome_order

        custom = ["chrA", "chrB", "chrC"]
        order = get_chromosome_order(species="canine", custom_order=custom)
        assert order == custom

    def test_unknown_species_raises_error(self):
        """Unknown species without custom order should raise."""
        from pylocuszoom.manhattan import get_chromosome_order

        with pytest.raises(ValueError, match="Unknown species"):
            get_chromosome_order(species="unknown")

    def test_no_species_or_custom_raises_error(self):
        """Must provide either species or custom_order."""
        from pylocuszoom.manhattan import get_chromosome_order

        with pytest.raises(ValueError, match="Must provide"):
            get_chromosome_order()

    def test_dog_alias_works(self):
        """'dog' should be an alias for 'canine'."""
        from pylocuszoom.manhattan import get_chromosome_order

        order = get_chromosome_order(species="dog")
        assert "1" in order
        assert "38" in order
        assert "X" in order

    def test_cat_alias_works(self):
        """'cat' should be an alias for 'feline'."""
        from pylocuszoom.manhattan import get_chromosome_order

        order = get_chromosome_order(species="cat")
        assert "A1" in order
        assert "E3" in order


class TestChromosomeColors:
    """Tests for chromosome color assignment."""

    def test_returns_correct_number_of_colors(self):
        """Should return exactly n colors."""
        from pylocuszoom.manhattan import get_chromosome_colors

        colors = get_chromosome_colors(10)
        assert len(colors) == 10

    def test_colors_are_hex_strings(self):
        """Colors should be hex strings."""
        from pylocuszoom.manhattan import get_chromosome_colors

        colors = get_chromosome_colors(5)
        for color in colors:
            assert color.startswith("#")
            assert len(color) == 7  # #RRGGBB format

    def test_colors_cycle_for_large_n(self):
        """Should handle more chromosomes than palette size."""
        from pylocuszoom.manhattan import get_chromosome_colors

        # Request more colors than typical palette size
        colors = get_chromosome_colors(100)
        assert len(colors) == 100


class TestPrepareManhattanData:
    """Tests for Manhattan data preparation."""

    @pytest.fixture
    def sample_gwas_df(self):
        """Create sample GWAS data for testing."""
        return pd.DataFrame(
            {
                "chrom": ["1", "1", "2", "2", "X"],
                "pos": [1000, 2000, 1500, 3000, 500],
                "p": [0.05, 1e-8, 0.01, 1e-10, 0.001],
                "snp": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            }
        )

    def test_adds_neg_log_p_column(self, sample_gwas_df):
        """Should compute -log10(p) for plotting."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        result = prepare_manhattan_data(sample_gwas_df, species="human")
        assert "_neg_log_p" in result.columns
        # Check calculation for first row
        expected = -np.log10(0.05)
        assert np.isclose(result["_neg_log_p"].iloc[0], expected, rtol=0.01)

    def test_adds_cumulative_position(self, sample_gwas_df):
        """Should compute cumulative x positions across chromosomes."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        result = prepare_manhattan_data(sample_gwas_df, species="human")
        assert "_cumulative_pos" in result.columns
        # Positions should be monotonically increasing when sorted
        sorted_result = result.sort_values("_cumulative_pos")
        assert sorted_result["_cumulative_pos"].is_monotonic_increasing

    def test_adds_color_column(self, sample_gwas_df):
        """Should assign colors to each variant."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        result = prepare_manhattan_data(sample_gwas_df, species="human")
        assert "_color" in result.columns
        assert all(c.startswith("#") for c in result["_color"])

    def test_stores_chrom_centers_in_attrs(self, sample_gwas_df):
        """Should store chromosome center positions for axis labels."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        result = prepare_manhattan_data(sample_gwas_df, species="human")
        assert "chrom_centers" in result.attrs
        centers = result.attrs["chrom_centers"]
        assert "1" in centers
        assert "2" in centers

    def test_validates_required_columns(self):
        """Should raise on missing required columns."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        df = pd.DataFrame({"wrong_col": [1, 2, 3]})
        with pytest.raises(ValueError, match="not found"):
            prepare_manhattan_data(df, species="human")

    def test_handles_integer_chromosomes(self):
        """Should handle int chromosome column."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        df = pd.DataFrame(
            {
                "chrom": [1, 1, 2],  # ints, not strings
                "pos": [1000, 2000, 1500],
                "p": [0.05, 0.01, 0.001],
            }
        )
        result = prepare_manhattan_data(df, species="human")
        assert len(result) == 3

    def test_handles_unknown_chromosomes(self):
        """Unknown chromosomes should be appended at end."""
        from pylocuszoom.manhattan import prepare_manhattan_data

        df = pd.DataFrame(
            {
                "chrom": ["1", "UNKNOWN"],
                "pos": [1000, 500],
                "p": [0.05, 0.01],
            }
        )
        result = prepare_manhattan_data(df, species="human")
        # Unknown chromosome should have larger cumulative position
        chrom1_pos = result[result["_chrom_str"] == "1"]["_cumulative_pos"].iloc[0]
        unknown_pos = result[result["_chrom_str"] == "UNKNOWN"]["_cumulative_pos"].iloc[
            0
        ]
        assert unknown_pos > chrom1_pos


class TestPrepareCategoricalData:
    """Tests for categorical (PheWAS-style) Manhattan data."""

    @pytest.fixture
    def sample_categorical_df(self):
        """Create sample categorical data."""
        return pd.DataFrame(
            {
                "phenotype": ["Height", "Weight", "BMI", "Height"],
                "p": [0.05, 1e-8, 0.01, 1e-10],
            }
        )

    def test_adds_neg_log_p(self, sample_categorical_df):
        """Should compute -log10(p)."""
        from pylocuszoom.manhattan import prepare_categorical_data

        result = prepare_categorical_data(
            sample_categorical_df, category_col="phenotype"
        )
        assert "_neg_log_p" in result.columns

    def test_adds_x_position(self, sample_categorical_df):
        """Should compute x positions for categories."""
        from pylocuszoom.manhattan import prepare_categorical_data

        result = prepare_categorical_data(
            sample_categorical_df, category_col="phenotype"
        )
        assert "_x_pos" in result.columns

    def test_adds_color(self, sample_categorical_df):
        """Should assign colors to categories."""
        from pylocuszoom.manhattan import prepare_categorical_data

        result = prepare_categorical_data(
            sample_categorical_df, category_col="phenotype"
        )
        assert "_color" in result.columns

    def test_validates_category_column(self):
        """Should raise on missing category column."""
        from pylocuszoom.manhattan import prepare_categorical_data

        df = pd.DataFrame({"p": [0.05]})
        with pytest.raises(ValueError, match="not found"):
            prepare_categorical_data(df, category_col="missing")

    def test_respects_custom_category_order(self, sample_categorical_df):
        """Custom order should be used."""
        from pylocuszoom.manhattan import prepare_categorical_data

        result = prepare_categorical_data(
            sample_categorical_df,
            category_col="phenotype",
            category_order=["BMI", "Weight", "Height"],
        )
        assert result.attrs["category_order"] == ["BMI", "Weight", "Height"]

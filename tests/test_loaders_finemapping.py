"""Tests for the fine-mapping result loaders."""

import pytest

from pylocuszoom.loaders import (
    load_caviar,
    load_finemap,
    load_polyfun,
    load_susie,
)
from pylocuszoom.schemas import LoaderValidationError


class TestSuSiELoader:
    """Tests for SuSiE file loader."""

    @pytest.fixture
    def susie_unset_cs_file(self, tmp_path):
        """Create a SuSiE file whose cs column holds -1 and a blank, never 0."""
        content = """pos\tpip\tcs\tsnp
1000000\t0.85\t1\trs123
1001000\t0.12\t-1\trs456
1002000\t0.02\t\trs789
1003000\t0.45\t2\trs101
"""
        filepath = tmp_path / "susie_unset.tsv"
        filepath.write_text(content)
        return filepath

    def test_load_susie_clamps_negative_credible_set(self, susie_unset_cs_file):
        """Test that a -1 credible set is clamped to 0."""
        df = load_susie(susie_unset_cs_file)

        assert df[df["pos"] == 1001000]["cs"].iloc[0] == 0  # Not -1, the file value

    def test_load_susie_fills_blank_credible_set(self, susie_unset_cs_file):
        """Test that a blank credible set cell becomes 0."""
        df = load_susie(susie_unset_cs_file)

        assert df[df["pos"] == 1002000]["cs"].iloc[0] == 0  # Not NaN, the file value

    def test_load_susie_basic(self, susie_file):
        """Test basic SuSiE file loading."""
        df = load_susie(susie_file)

        assert "pos" in df.columns
        assert "pip" in df.columns
        assert "cs" in df.columns
        assert len(df) == 4

    def test_load_susie_credible_sets(self, susie_file):
        """Test credible set values."""
        df = load_susie(susie_file)

        # Check credible set assignments
        assert df[df["pos"] == 1000000]["cs"].iloc[0] == 1
        assert df[df["pos"] == 1002000]["cs"].iloc[0] == 0  # Not in CS
        assert df[df["pos"] == 1003000]["cs"].iloc[0] == 2


class TestFINEMAPLoader:
    """Tests for FINEMAP file loader."""

    @pytest.fixture
    def finemap_file(self, tmp_path):
        """Create a temporary FINEMAP .snp output file."""
        content = """index rsid chromosome position allele1 allele2 maf beta se z prob log10bf
1 rs123 1 1000000 A G 0.3 0.5 0.2 2.5 0.85 2.1
2 rs456 1 1001000 C T 0.2 0.3 0.15 2.0 0.12 1.5
3 rs789 1 1002000 G A 0.4 -0.2 0.1 -2.0 0.02 0.5
4 rs101 1 1003000 T C 0.1 0.4 0.25 1.6 0.01 0.3
"""
        filepath = tmp_path / "results.snp"
        filepath.write_text(content)
        return filepath

    def test_load_finemap_basic(self, finemap_file):
        """Test basic FINEMAP file loading."""
        df = load_finemap(finemap_file)

        assert "pos" in df.columns
        assert "pip" in df.columns
        assert len(df) == 4

    def test_load_finemap_assigns_credible_set(self, finemap_file):
        """Test that FINEMAP loader assigns credible sets based on cumsum."""
        df = load_finemap(finemap_file)

        # Sorted by PIP: 0.85 + 0.12 = 0.97 > 0.95, so first 2 in CS
        assert "cs" in df.columns
        # The first variant in sorted order (0.85) should be in credible set
        cs_variants = df[df["cs"] == 1]
        assert len(cs_variants) >= 1

    def test_load_finemap_values(self, finemap_file):
        """Test that values are loaded correctly."""
        df = load_finemap(finemap_file)

        # After loading, find the variant with rsid rs123
        rs123 = df[df["rs"] == "rs123"]
        if len(rs123) > 0:
            assert rs123["pip"].iloc[0] == 0.85


class TestCAVIARLoader:
    """Tests for CAVIAR file loader."""

    @pytest.fixture
    def caviar_file(self, tmp_path):
        """Create a temporary CAVIAR .set output file."""
        content = """rs123 0.85
rs456 0.12
rs789 0.02
rs101 0.01
"""
        filepath = tmp_path / "results.set"
        filepath.write_text(content)
        return filepath

    def test_load_caviar_basic(self, caviar_file):
        """Test basic CAVIAR file loading."""
        df = load_caviar(caviar_file)

        assert "rs" in df.columns
        assert "pip" in df.columns
        assert len(df) == 4

    def test_load_caviar_assigns_credible_set(self, caviar_file):
        """Test that CAVIAR loader assigns credible sets."""
        df = load_caviar(caviar_file)

        assert "cs" in df.columns
        # Top variants should be in credible set (cumsum <= 0.95)
        assert df[df["rs"] == "rs123"]["cs"].iloc[0] == 1

    def test_load_caviar_no_position_column(self, caviar_file):
        """Test that CAVIAR output doesn't include position column."""
        df = load_caviar(caviar_file)

        # CAVIAR doesn't include positions - user needs to merge
        assert "pos" not in df.columns
        # But should have rs and pip
        assert "rs" in df.columns
        assert "pip" in df.columns


class TestPolyFunLoader:
    """Tests for PolyFun file loader."""

    @pytest.fixture
    def polyfun_file(self, tmp_path):
        """Create a temporary PolyFun output file."""
        content = """CHR BP SNP A1 A2 PIP CREDIBLE_SET BETA SE
1 1000000 rs123 A G 0.85 1 0.5 0.2
1 1001000 rs456 C T 0.12 1 0.3 0.15
1 1002000 rs789 G A 0.02 0 -0.2 0.1
"""
        filepath = tmp_path / "polyfun.txt"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def polyfun_file_no_pip(self, tmp_path):
        """Create a PolyFun file with no column the spec can map to pip."""
        content = """CHR BP SNP A1 A2 CREDIBLE_SET BETA SE
1 1000000 rs123 A G 1 0.5 0.2
1 1001000 rs456 C T 1 0.3 0.15
1 1002000 rs789 G A 0 -0.2 0.1
"""
        filepath = tmp_path / "polyfun_no_pip.txt"
        filepath.write_text(content)
        return filepath

    def test_load_polyfun_basic(self, polyfun_file):
        """Test basic PolyFun file loading."""
        df = load_polyfun(polyfun_file)

        assert "pos" in df.columns
        assert "pip" in df.columns
        assert "cs" in df.columns
        assert len(df) == 3

    def test_load_polyfun_preserves_credible_set(self, polyfun_file):
        """Test that PolyFun loader preserves CREDIBLE_SET column."""
        df = load_polyfun(polyfun_file)

        assert df[df["pos"] == 1000000]["cs"].iloc[0] == 1
        assert df[df["pos"] == 1002000]["cs"].iloc[0] == 0

    def test_load_polyfun_unmappable_pip_raises(self, polyfun_file_no_pip):
        """An unmappable pip column fails at load time, naming the column."""
        with pytest.raises(LoaderValidationError) as exc_info:
            load_polyfun(polyfun_file_no_pip)

        assert "pip" in str(exc_info.value)

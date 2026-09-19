"""Tests for the eQTL file loaders."""

import pytest

from pylocuszoom.loaders import (
    load_eqtl_catalogue,
    load_gtex_eqtl,
    load_matrixeqtl,
)


class TestGTExEQTLLoader:
    """Tests for GTEx eQTL file loader."""

    @pytest.fixture
    def gtex_file(self, tmp_path):
        """Create a temporary GTEx eQTL file."""
        content = """variant_id\tgene_id\tpval_nominal\tslope
chr1_1000000_A_G_b38\tENSG00001\t1e-6\t0.5
chr1_1001000_C_T_b38\tENSG00001\t0.01\t-0.3
"""
        filepath = tmp_path / "gtex.txt"
        filepath.write_text(content)
        return filepath

    def test_load_gtex_effect_size_column(self, gtex_file):
        """Test that GTEx loader outputs effect_size column (not effect)."""
        df = load_gtex_eqtl(gtex_file)

        # Should standardize to effect_size for compatibility with plotter
        assert "effect_size" in df.columns
        assert df["effect_size"].iloc[0] == 0.5


class TestEQTLCatalogueLoader:
    """Tests for eQTL Catalogue file loader."""

    @pytest.fixture
    def eqtl_catalogue_file(self, tmp_path):
        """Create a temporary eQTL Catalogue file."""
        content = """molecular_trait_id\tgene_id\tvariant\tchromosome\tposition\tref\talt\tac\tan\tmaf\tbeta\tse\tpvalue
ENSG00001_ENST00001\tENSG00001\t1_1000000_A_G\t1\t1000000\tA\tG\t100\t1000\t0.1\t0.5\t0.1\t1e-6
ENSG00001_ENST00001\tENSG00001\t1_1001000_C_T\t1\t1001000\tC\tT\t200\t1000\t0.2\t-0.3\t0.15\t0.01
"""
        filepath = tmp_path / "eqtl_catalogue.tsv"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def eqtl_catalogue_substring_file(self, tmp_path):
        """Create an eQTL Catalogue file whose genes TP5 and TP53 overlap."""
        content = """chromosome\tposition\tgene_id\tbeta\tpvalue
1\t1000000\tTP5\t0.5\t1e-6
1\t1001000\tTP53\t-0.3\t0.01
1\t1002000\tBRCA1\t0.2\t0.02
"""
        filepath = tmp_path / "eqtl_catalogue_substring.tsv"
        filepath.write_text(content)
        return filepath

    def test_load_eqtl_catalogue_basic(self, eqtl_catalogue_file):
        """Test basic eQTL Catalogue file loading."""
        df = load_eqtl_catalogue(eqtl_catalogue_file)

        assert "pos" in df.columns
        assert "p_value" in df.columns
        assert "gene" in df.columns
        assert len(df) == 2

    def test_load_eqtl_catalogue_effect_size_column(self, eqtl_catalogue_file):
        """Test that eQTL Catalogue loader outputs effect_size column."""
        df = load_eqtl_catalogue(eqtl_catalogue_file)

        assert "effect_size" in df.columns
        assert df["effect_size"].iloc[0] == 0.5

    def test_load_eqtl_catalogue_gene_filter(self, eqtl_catalogue_file):
        """Test gene filtering works."""
        df = load_eqtl_catalogue(eqtl_catalogue_file, gene="ENSG00001")

        assert len(df) == 2

    def test_load_eqtl_catalogue_gene_filter_is_substring(
        self, eqtl_catalogue_substring_file
    ):
        """Test that eQTL Catalogue gene filtering matches on substring."""
        df = load_eqtl_catalogue(eqtl_catalogue_substring_file, gene="TP5")

        assert len(df) == 2  # Not 1; "contains" also matches TP53
        assert set(df["gene"]) == {"TP5", "TP53"}


class TestMatrixEQTLLoader:
    """Tests for MatrixEQTL file loader."""

    @pytest.fixture
    def matrixeqtl_file(self, tmp_path):
        """Create a temporary MatrixEQTL output file."""
        content = """SNP\tgene\tbeta\tt-stat\tp-value\tFDR
rs123\tBRCA1\t0.5\t3.5\t1e-6\t1e-5
rs456\tBRCA1\t-0.3\t-2.1\t0.03\t0.1
rs789\tTP53\t0.2\t2.0\t0.04\t0.12
"""
        filepath = tmp_path / "matrixeqtl.txt"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def matrixeqtl_substring_file(self, tmp_path):
        """Create a MatrixEQTL file whose genes TP5 and TP53 overlap."""
        content = """SNP\tgene\tbeta\tt-stat\tp-value\tFDR
rs123\tTP5\t0.5\t3.5\t1e-6\t1e-5
rs456\tTP53\t-0.3\t-2.1\t0.03\t0.1
rs789\tBRCA1\t0.2\t2.0\t0.04\t0.12
"""
        filepath = tmp_path / "matrixeqtl_substring.txt"
        filepath.write_text(content)
        return filepath

    def test_load_matrixeqtl_basic(self, matrixeqtl_file):
        """Test basic MatrixEQTL file loading."""
        df = load_matrixeqtl(matrixeqtl_file)

        assert "rs" in df.columns
        assert "gene" in df.columns
        assert "p_value" in df.columns
        assert len(df) == 3

    def test_load_matrixeqtl_effect_size_column(self, matrixeqtl_file):
        """Test that MatrixEQTL loader outputs effect_size column."""
        df = load_matrixeqtl(matrixeqtl_file)

        assert "effect_size" in df.columns
        assert df["effect_size"].iloc[0] == 0.5

    def test_load_matrixeqtl_gene_filter(self, matrixeqtl_file):
        """Test gene filtering works."""
        df = load_matrixeqtl(matrixeqtl_file, gene="BRCA1")

        assert len(df) == 2
        assert all(df["gene"] == "BRCA1")

    def test_load_matrixeqtl_gene_filter_is_exact(self, matrixeqtl_substring_file):
        """Test that MatrixEQTL gene filtering matches on equality."""
        df = load_matrixeqtl(matrixeqtl_substring_file, gene="TP5")

        assert len(df) == 1  # Not 2; "exact" excludes TP53
        assert set(df["gene"]) == {"TP5"}


class TestGTExCoordinates:
    def test_loaded_chromosomes_keep_equal_positions_distinct(self, tmp_path):
        import pandas as pd

        from pylocuszoom.eqtl import (
            calculate_colocalization_overlap,
            filter_eqtl_by_region,
        )

        path = tmp_path / "gtex.tsv"
        variants = ["chr1_1000_A_G_b38", "chr2_1000_A_G_b38"]
        pd.DataFrame(
            {
                "variant_id": variants,
                "gene_id": ["GENE1", "GENE1"],
                "pval_nominal": [1e-8, 1e-9],
            }
        ).to_csv(path, sep="\t", index=False)
        loaded = load_gtex_eqtl(path)
        assert loaded["chr"].tolist() == ["1", "2"]
        assert loaded["variant_id"].tolist() == variants
        region = filter_eqtl_by_region(loaded, 1, 900, 1100)
        assert region["variant_id"].tolist() == [variants[0]]
        gwas = pd.DataFrame({"chr": [2], "pos": [1000], "p_value": [1e-8]})
        overlap = calculate_colocalization_overlap(gwas, region)
        assert overlap.empty

    def test_relative_tss_distance_is_not_an_absolute_position(self, tmp_path):
        from pylocuszoom.exceptions import LoaderValidationError

        path = tmp_path / "relative.tsv"
        path.write_text("tss_distance\tpval_nominal\tgene_id\n500\t1e-8\tGENE1\n")
        with pytest.raises(LoaderValidationError, match="pos"):
            load_gtex_eqtl(path)

    @pytest.mark.parametrize("variant", ["bad", "chr1_bad_A_G_b38", None])
    def test_invalid_variant_id_raises_loader_error(self, tmp_path, variant):
        import pandas as pd

        from pylocuszoom.exceptions import LoaderValidationError

        path = tmp_path / "invalid.tsv"
        pd.DataFrame({"variant_id": [variant], "pval_nominal": [1e-8]}).to_csv(
            path, sep="\t", index=False
        )
        with pytest.raises(LoaderValidationError, match="variant_id"):
            load_gtex_eqtl(path)

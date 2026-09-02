"""Tests for file format loaders and pydantic validation."""

import io

import pandas as pd
import pytest

from pylocuszoom.loaders import (
    load_bed,
    load_bolt_lmm,
    load_caviar,
    load_ensembl_genes,
    load_eqtl_catalogue,
    load_finemap,
    load_gemma,
    load_gtex_eqtl,
    load_gtf,
    load_gwas,
    load_gwas_catalog,
    load_matrixeqtl,
    load_plink_assoc,
    load_polyfun,
    load_regenie,
    load_saige,
    load_susie,
)
from pylocuszoom.schemas import (
    LoaderValidationError,
    validate_eqtl_dataframe,
    validate_finemapping_dataframe,
    validate_genes_dataframe,
    validate_gwas_dataframe,
)

# =============================================================================
# Fixtures for test data files
# =============================================================================


@pytest.fixture
def plink_assoc_file(tmp_path):
    """Create a temporary PLINK .assoc file."""
    content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
1 rs456 1001000 G ADD 1000 0.3 1.5 0.1
1 rs789 1002000 T ADD 1000 -0.2 -1.0 1e-8
"""
    filepath = tmp_path / "test.assoc"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def regenie_file(tmp_path):
    """Create a temporary REGENIE file."""
    content = """CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA
1 1000000 rs123 A G 0.3 1000 ADD 0.5 0.2 6.25 2.0 NA
1 1001000 rs456 C T 0.2 1000 ADD 0.3 0.15 4.0 1.5 NA
1 1002000 rs789 G A 0.4 1000 ADD -0.2 0.1 4.0 8.0 NA
"""
    filepath = tmp_path / "test.regenie"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def susie_file(tmp_path):
    """Create a temporary SuSiE results file."""
    content = """pos\tpip\tcs\tsnp
1000000\t0.85\t1\trs123
1001000\t0.12\t1\trs456
1002000\t0.02\t0\trs789
1003000\t0.45\t2\trs101
"""
    filepath = tmp_path / "susie.tsv"
    filepath.write_text(content)
    return filepath


@pytest.fixture
def bed_file(tmp_path):
    """Create a temporary BED file."""
    content = """chr1\t1000000\t1020000\tGENE1
chr1\t1050000\t1080000\tGENE2
chr1\t1100000\t1150000\tGENE3
"""
    filepath = tmp_path / "genes.bed"
    filepath.write_text(content)
    return filepath


# =============================================================================
# Test GWAS Loaders
# =============================================================================


class TestPLINKLoader:
    """Tests for PLINK association file loader."""

    @pytest.fixture
    def plink_glm_hash_file(self, tmp_path):
        """Create a PLINK2 --glm file with its native '#CHROM' header line."""
        content = """#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP
1\t1000000\trs123\tA\tG\tG\tADD\t1000\t0.5\t0.2\t2.5\t0.01
1\t1001000\trs456\tC\tT\tT\tADD\t1000\t0.3\t0.15\t1.5\t0.1
1\t1002000\trs789\tG\tA\tA\tADD\t1000\t-0.2\t0.1\t-1.0\t1e-8
"""
        filepath = tmp_path / "glm_hash.assoc.linear"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def plink_glm_file(self, tmp_path):
        """Create a PLINK2 --glm file with the leading '#' stripped."""
        content = """CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP
1\t1000000\trs123\tA\tG\tG\tADD\t1000\t0.5\t0.2\t2.5\t0.01
1\t1001000\trs456\tC\tT\tT\tADD\t1000\t0.3\t0.15\t1.5\t0.1
1\t1002000\trs789\tG\tA\tA\tADD\t1000\t-0.2\t0.1\t-1.0\t1e-8
"""
        filepath = tmp_path / "glm.assoc.linear"
        filepath.write_text(content)
        return filepath

    def test_load_plink_glm_native_hash_header(self, plink_glm_hash_file):
        """Load a PLINK2 --glm file with its native '#CHROM' header line.

        pandas treats '#' as a comment prefix, so a spec-level comment="#"
        swallowed the header and promoted the first variant to column labels.
        The chr candidate '#CHROM' exists precisely for this header.
        """
        df = load_plink_assoc(plink_glm_hash_file)

        assert len(df) == 3, "the header line must not be counted as a variant"
        assert df["chr"].iloc[0] == 1
        assert df["ps"].iloc[0] == 1000000
        assert df["rs"].iloc[0] == "rs123"
        assert df["p_wald"].iloc[0] == 0.01

    def test_load_gwas_detects_plink2_glm(self, plink_glm_hash_file):
        """A .glm.linear name routes to the PLINK loader, not the warn-default."""
        renamed = plink_glm_hash_file.parent / "plink2.PHENO1.glm.linear"
        plink_glm_hash_file.rename(renamed)

        df = load_gwas(renamed)

        assert len(df) == 3
        assert df["rs"].iloc[0] == "rs123"

    def test_load_plink_glm_style_aliases_without_hash(self, plink_glm_file):
        """Map PLINK2 --glm aliases when no '#' hides the header row."""
        df = load_plink_assoc(plink_glm_file)

        assert len(df) == 3  # Not 2; no '#' line for comment="#" to swallow
        assert df["ps"].iloc[0] == 1000000  # POS, the second pos_col candidate
        assert df["rs"].iloc[0] == "rs123"  # ID, the second rs_col candidate
        assert df["chr"].iloc[0] == 1  # CHROM, the third chr candidate
        assert df["p_wald"].iloc[0] == 0.01

    def test_load_plink_assoc_basic(self, plink_assoc_file):
        """Test basic PLINK file loading."""
        df = load_plink_assoc(plink_assoc_file)

        assert "ps" in df.columns
        assert "p_wald" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_plink_assoc_custom_columns(self, plink_assoc_file):
        """Test PLINK loader with custom column names."""
        df = load_plink_assoc(
            plink_assoc_file,
            pos_col="position",
            p_col="pvalue",
            rs_col="snp_id",
        )

        assert "position" in df.columns
        assert "pvalue" in df.columns
        assert "snp_id" in df.columns

    def test_load_plink_assoc_values_correct(self, plink_assoc_file):
        """Test that loaded values are correct."""
        df = load_plink_assoc(plink_assoc_file)

        assert df["ps"].iloc[0] == 1000000
        assert df["p_wald"].iloc[0] == 0.01
        assert df["rs"].iloc[0] == "rs123"


class TestREGENIELoader:
    """Tests for REGENIE file loader."""

    @pytest.fixture
    def regenie_p_only_file(self, tmp_path):
        """Create a REGENIE file carrying a raw P column and no LOG10P."""
        content = """CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ P
1 1000000 rs123 A G 0.3 1000 ADD 0.5 0.2 6.25 0.004
1 1001000 rs456 C T 0.2 1000 ADD 0.3 0.15 4.0 0.05
1 1002000 rs789 G A 0.4 1000 ADD -0.2 0.1 4.0 2e-9
"""
        filepath = tmp_path / "test_p_only.regenie"
        filepath.write_text(content)
        return filepath

    def test_load_regenie_renames_raw_p_column(self, regenie_p_only_file):
        """Test that REGENIE renames P when LOG10P is absent."""
        df = load_regenie(regenie_p_only_file)

        assert df["p_wald"].iloc[0] == 0.004  # Taken as-is, not 10**(-0.004)
        assert df["p_wald"].iloc[2] == pytest.approx(2e-9, rel=0.01)

    def test_load_regenie_basic(self, regenie_file):
        """Test basic REGENIE file loading."""
        df = load_regenie(regenie_file)

        assert "ps" in df.columns
        assert "p_wald" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_regenie_log10p_conversion(self, regenie_file):
        """Test that LOG10P is converted to p-value."""
        df = load_regenie(regenie_file)

        # LOG10P=2.0 -> p=0.01, LOG10P=8.0 -> p=1e-8
        assert df["p_wald"].iloc[0] == pytest.approx(0.01, rel=0.01)
        assert df["p_wald"].iloc[2] == pytest.approx(1e-8, rel=0.01)


class TestSAIGELoader:
    """Tests for SAIGE file loader."""

    @pytest.fixture
    def saige_file(self, tmp_path):
        """Create a temporary SAIGE file with both p.value columns."""
        content = """CHR\tPOS\tMarkerID\tp.value\tp.value.NA
1\t1000000\trs123\t0.01\t0.005
1\t1001000\trs456\t0.05\t0.04
1\t1002000\trs789\t1e-8\t1e-9
"""
        filepath = tmp_path / "test.saige"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def saige_file_no_spa(self, tmp_path):
        """Create a SAIGE file with only p.value (no SPA-adjusted)."""
        content = """CHR\tPOS\tMarkerID\tp.value
1\t1000000\trs123\t0.01
1\t1001000\trs456\t0.05
"""
        filepath = tmp_path / "test_nospa.saige"
        filepath.write_text(content)
        return filepath

    def test_load_saige_prefers_spa_adjusted(self, saige_file):
        """Test that SAIGE loader prefers p.value.NA (SPA-adjusted) over p.value."""
        df = load_saige(saige_file)

        assert "p_wald" in df.columns
        # Should use SPA-adjusted p-values (p.value.NA column)
        assert df["p_wald"].iloc[0] == 0.005  # Not 0.01
        assert df["p_wald"].iloc[2] == pytest.approx(1e-9, rel=0.01)

    def test_load_saige_fallback_to_pvalue(self, saige_file_no_spa):
        """Test that SAIGE loader falls back to p.value when p.value.NA missing."""
        df = load_saige(saige_file_no_spa)

        assert "p_wald" in df.columns
        assert df["p_wald"].iloc[0] == 0.01

    def test_load_saige_basic(self, saige_file):
        """Test basic SAIGE file loading."""
        df = load_saige(saige_file)

        assert "ps" in df.columns
        assert "p_wald" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3


class TestBOLTLMMLoader:
    """Tests for BOLT-LMM file loader."""

    @pytest.fixture
    def bolt_file(self, tmp_path):
        """Create a temporary BOLT-LMM stats file."""
        content = """SNP\tCHR\tBP\tGENPOS\tALLELE1\tALLELE0\tA1FREQ\tF_MISS\tBETA\tSE\tP_BOLT_LMM_INF\tP_BOLT_LMM
rs123\t1\t1000000\t0.01\tA\tG\t0.3\t0.01\t0.5\t0.2\t0.01\t0.005
rs456\t1\t1001000\t0.02\tC\tT\t0.2\t0.02\t0.3\t0.15\t0.05\t0.04
rs789\t1\t1002000\t0.03\tG\tA\t0.4\t0.01\t-0.2\t0.1\t1e-8\t1e-9
"""
        filepath = tmp_path / "test.stats"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def bolt_file_inf_only(self, tmp_path):
        """Create BOLT-LMM file with only infinitesimal p-value."""
        content = """SNP\tCHR\tBP\tALLELE1\tALLELE0\tA1FREQ\tBETA\tSE\tP_BOLT_LMM_INF
rs123\t1\t1000000\tA\tG\t0.3\t0.5\t0.2\t0.01
"""
        filepath = tmp_path / "test_inf.stats"
        filepath.write_text(content)
        return filepath

    def test_load_bolt_basic(self, bolt_file):
        """Test basic BOLT-LMM file loading."""
        df = load_bolt_lmm(bolt_file)

        assert "ps" in df.columns
        assert "p_wald" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_bolt_prefers_full_model(self, bolt_file):
        """Test that BOLT-LMM loader prefers P_BOLT_LMM over P_BOLT_LMM_INF."""
        df = load_bolt_lmm(bolt_file)

        # Should use P_BOLT_LMM (full model), not P_BOLT_LMM_INF
        assert df["p_wald"].iloc[0] == 0.005  # Not 0.01
        assert df["p_wald"].iloc[2] == pytest.approx(1e-9, rel=0.01)

    def test_load_bolt_fallback_to_inf(self, bolt_file_inf_only):
        """Test that BOLT-LMM loader falls back to P_BOLT_LMM_INF."""
        df = load_bolt_lmm(bolt_file_inf_only)

        assert "p_wald" in df.columns
        assert df["p_wald"].iloc[0] == 0.01


class TestGEMMALoader:
    """Tests for GEMMA file loader."""

    @pytest.fixture
    def gemma_file(self, tmp_path):
        """Create a temporary GEMMA .assoc.txt file."""
        content = """chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tlogl_H1\tl_remle\tp_wald
1\trs123\t1000000\t0\tA\tG\t0.3\t0.5\t0.2\t100\t0.5\t0.01
1\trs456\t1001000\t0\tC\tT\t0.2\t0.3\t0.15\t95\t0.4\t0.05
1\trs789\t1002000\t0\tG\tA\t0.4\t-0.2\t0.1\t110\t0.6\t1e-8
"""
        filepath = tmp_path / "output.assoc.txt"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def gemma_lrt_file(self, tmp_path):
        """Create GEMMA file with p_lrt instead of p_wald."""
        content = """chr\trs\tps\tallele1\tallele0\taf\tbeta\tse\tp_lrt
1\trs123\t1000000\tA\tG\t0.3\t0.5\t0.2\t0.02
"""
        filepath = tmp_path / "output_lrt.assoc.txt"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def gemma_all_p_file(self, tmp_path):
        """Create a GEMMA file carrying p_wald, p_lrt and p_score together."""
        content = """chr\trs\tps\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score
1\trs123\t1000000\tA\tG\t0.3\t0.5\t0.2\t0.001\t0.02\t0.3
1\trs456\t1001000\tC\tT\t0.2\t0.3\t0.15\t0.004\t0.06\t0.4
"""
        filepath = tmp_path / "output_all_p.assoc.txt"
        filepath.write_text(content)
        return filepath

    def test_load_gemma_prefers_wald_over_lrt_and_score(self, gemma_all_p_file):
        """Test that GEMMA prefers p_wald over p_lrt and p_score."""
        df = load_gemma(gemma_all_p_file)

        assert df["p_wald"].iloc[0] == 0.001  # Not 0.02 (p_lrt), not 0.3 (p_score)

    def test_load_gemma_precedence_with_custom_p_col(self, gemma_all_p_file):
        """Test GEMMA p-value precedence through a custom output column.

        With the default p_col a reordered candidate tuple would rename p_lrt
        onto the existing p_wald column, so the assertion would hit duplicate
        labels rather than a clean value mismatch.
        """
        df = load_gemma(gemma_all_p_file, p_col="pval")

        assert df["pval"].iloc[0] == 0.001  # Not 0.02 (p_lrt), not 0.3 (p_score)

    def test_load_gemma_basic(self, gemma_file):
        """Test basic GEMMA file loading."""
        df = load_gemma(gemma_file)

        assert "ps" in df.columns
        assert "p_wald" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_gemma_values(self, gemma_file):
        """Test that values are loaded correctly."""
        df = load_gemma(gemma_file)

        assert df["ps"].iloc[0] == 1000000
        assert df["p_wald"].iloc[0] == 0.01

    def test_load_gemma_fallback_to_lrt(self, gemma_lrt_file):
        """Test GEMMA loader falls back to p_lrt when p_wald missing."""
        df = load_gemma(gemma_lrt_file)

        assert "p_wald" in df.columns
        assert df["p_wald"].iloc[0] == 0.02


class TestGWASCatalogLoader:
    """Tests for GWAS Catalog file loader."""

    @pytest.fixture
    def catalog_file(self, tmp_path):
        """Create a temporary GWAS Catalog file."""
        content = """chromosome\tbase_pair_location\tvariant_id\tp_value\tbeta
1\t1000000\trs123\t0.01\t0.5
1\t1001000\trs456\t0.001\t0.3
1\t1002000\trs789\t1e-8\t-0.2
"""
        filepath = tmp_path / "gwas_catalog.tsv"
        filepath.write_text(content)
        return filepath

    def test_load_gwas_catalog_basic(self, catalog_file):
        """Test basic GWAS Catalog file loading."""
        df = load_gwas_catalog(catalog_file)

        assert "ps" in df.columns
        assert "p_wald" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_gwas_catalog_values(self, catalog_file):
        """Test that values are mapped correctly."""
        df = load_gwas_catalog(catalog_file)

        assert df["ps"].iloc[0] == 1000000
        assert df["p_wald"].iloc[0] == 0.01
        assert df["rs"].iloc[0] == "rs123"


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


class TestAutoFormatDetection:
    """Tests for automatic format detection."""

    def test_load_gwas_detects_plink(self, plink_assoc_file):
        """Test that load_gwas auto-detects PLINK format."""
        df = load_gwas(plink_assoc_file)
        assert "ps" in df.columns
        assert len(df) == 3

    def test_load_gwas_detects_regenie(self, regenie_file):
        """Test that load_gwas auto-detects REGENIE format."""
        df = load_gwas(regenie_file)
        assert "ps" in df.columns
        assert len(df) == 3

    def test_load_gwas_detects_bolt_from_stats_ext(self, tmp_path):
        """Test that .stats extension is detected as BOLT-LMM format."""
        content = """SNP\tCHR\tBP\tALLELE1\tALLELE0\tA1FREQ\tBETA\tSE\tP_BOLT_LMM_INF
rs123\t1\t1000000\tA\tG\t0.3\t0.5\t0.2\t0.01
"""
        filepath = tmp_path / "result.stats"
        filepath.write_text(content)
        df = load_gwas(filepath)
        assert "ps" in df.columns
        assert len(df) == 1

    def test_load_gwas_unknown_format_warns_and_defaults_plink(self, tmp_path):
        """Unknown file extension logs warning and defaults to plink loader."""
        content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
"""
        filepath = tmp_path / "results.unknown_ext"
        filepath.write_text(content)

        from pylocuszoom.logging import enable_logging

        log_capture = io.StringIO()
        enable_logging("WARNING", sink=log_capture)
        try:
            df = load_gwas(filepath)
        finally:
            enable_logging("INFO")

        assert "ps" in df.columns
        assert len(df) == 1
        message = log_capture.getvalue()
        assert "results.unknown_ext" in message
        assert "%s" not in message

    def test_load_gwas_explicit_format_overrides_detection(self, tmp_path):
        """Explicit format= parameter overrides auto-detection."""
        # File with .stats extension but content is PLINK format
        content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
"""
        filepath = tmp_path / "tricky.stats"
        filepath.write_text(content)
        # Force plink format despite .stats extension
        df = load_gwas(filepath, format="plink")
        assert "ps" in df.columns

    def test_load_gwas_invalid_format_raises(self):
        """Invalid format name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown format"):
            load_gwas("/fake/path", format="nonexistent_format")

    def test_load_gwas_detects_assoc_linear(self, tmp_path):
        """Test .assoc.linear extension is detected as plink."""
        content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
"""
        filepath = tmp_path / "result.assoc.linear"
        filepath.write_text(content)
        df = load_gwas(filepath)
        assert "ps" in df.columns

    def test_load_gwas_detects_qassoc(self, tmp_path):
        """Test .qassoc extension is detected as plink."""
        content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
"""
        filepath = tmp_path / "result.qassoc"
        filepath.write_text(content)
        df = load_gwas(filepath)
        assert "ps" in df.columns


# =============================================================================
# Test Fine-mapping Loaders
# =============================================================================


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


# =============================================================================
# Test Gene Annotation Loaders
# =============================================================================


class TestBEDLoader:
    """Tests for BED file loader."""

    def test_load_bed_basic(self, bed_file):
        """Test basic BED file loading."""
        df = load_bed(bed_file)

        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "gene_name" in df.columns
        assert len(df) == 3

    def test_load_bed_chromosome_cleaned(self, bed_file):
        """Test that chromosome prefix is removed."""
        df = load_bed(bed_file)

        # "chr1" should become "1"
        assert df["chr"].iloc[0] == "1"

    def test_load_bed12_format(self, tmp_path):
        """Test BED12 format with extra columns doesn't break."""
        # BED12 has: chr, start, end, name, score, strand, thickStart, thickEnd,
        #            itemRgb, blockCount, blockSizes, blockStarts
        content = """chr1\t1000000\t1020000\tGENE1\t100\t+\t1000500\t1019500\t0\t3\t100,200,300\t0,5000,19700
chr1\t1050000\t1080000\tGENE2\t200\t-\t1050000\t1080000\t0\t2\t500,600\t0,29400
"""
        filepath = tmp_path / "genes.bed12"
        filepath.write_text(content)

        df = load_bed(filepath)

        assert len(df) == 2
        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "gene_name" in df.columns
        assert df["gene_name"].iloc[0] == "GENE1"

    def test_load_bed_7_columns(self, tmp_path):
        """Test BED with 7 columns (more than 6)."""
        content = """chr1\t1000000\t1020000\tGENE1\t100\t+\textra_col
chr1\t1050000\t1080000\tGENE2\t200\t-\tmore_data
"""
        filepath = tmp_path / "genes.bed7"
        filepath.write_text(content)

        df = load_bed(filepath)

        assert len(df) == 2
        assert "gene_name" in df.columns
        assert df["gene_name"].iloc[0] == "GENE1"


class TestGTFLoader:
    """Tests for GTF file loader."""

    @pytest.fixture
    def gtf_file(self, tmp_path):
        """Create a temporary GTF file with gene_name before gene_id."""
        # Note: loader extracts first matching attribute (gene_name or gene_id)
        # so we put gene_name first to test gene_name extraction
        content = """##description: test GTF file
chr1\tENSEMBL\tgene\t1000000\t1020000\t.\t+\t.\tgene_name "BRCA1"; gene_id "ENSG00001"; gene_biotype "protein_coding";
chr1\tENSEMBL\texon\t1000000\t1005000\t.\t+\t.\tgene_name "BRCA1"; gene_id "ENSG00001"; exon_number 1;
chr1\tENSEMBL\texon\t1015000\t1020000\t.\t+\t.\tgene_name "BRCA1"; gene_id "ENSG00001"; exon_number 2;
chr1\tENSEMBL\tgene\t1050000\t1080000\t.\t-\t.\tgene_name "TP53"; gene_id "ENSG00002"; gene_biotype "protein_coding";
"""
        filepath = tmp_path / "genes.gtf"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def gtf_file_gene_id_only(self, tmp_path):
        """Create a GTF file with only gene_id (no gene_name)."""
        content = """##description: test GTF file
chr1\tENSEMBL\tgene\t1000000\t1020000\t.\t+\t.\tgene_id "ENSG00001"; gene_biotype "protein_coding";
"""
        filepath = tmp_path / "genes_id_only.gtf"
        filepath.write_text(content)
        return filepath

    def test_load_gtf_genes(self, gtf_file):
        """Test loading genes from GTF file."""
        df = load_gtf(gtf_file, feature_type="gene")

        assert len(df) == 2
        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "gene_name" in df.columns
        assert "strand" in df.columns

    def test_load_gtf_exons(self, gtf_file):
        """Test loading exons from GTF file."""
        df = load_gtf(gtf_file, feature_type="exon")

        assert len(df) == 2
        assert all(df["gene_name"] == "BRCA1")

    def test_load_gtf_chromosome_cleaned(self, gtf_file):
        """Test that chromosome prefix is removed."""
        df = load_gtf(gtf_file, feature_type="gene")

        # "chr1" should become "1"
        assert df["chr"].iloc[0] == "1"

    def test_load_gtf_extracts_gene_name(self, gtf_file):
        """Test that gene_name is extracted from attributes."""
        df = load_gtf(gtf_file, feature_type="gene")

        assert "BRCA1" in df["gene_name"].values
        assert "TP53" in df["gene_name"].values

    def test_load_gtf_preserves_strand(self, gtf_file):
        """Test that strand information is preserved."""
        df = load_gtf(gtf_file, feature_type="gene")

        brca1 = df[df["gene_name"] == "BRCA1"]
        tp53 = df[df["gene_name"] == "TP53"]
        assert brca1["strand"].iloc[0] == "+"
        assert tp53["strand"].iloc[0] == "-"

    def test_load_gtf_fallback_to_gene_id(self, gtf_file_gene_id_only):
        """Test that loader falls back to gene_id when gene_name missing."""
        df = load_gtf(gtf_file_gene_id_only, feature_type="gene")

        assert len(df) == 1
        assert df["gene_name"].iloc[0] == "ENSG00001"


class TestEnsemblGenesLoader:
    """Tests for Ensembl BioMart gene export loader."""

    @pytest.fixture
    def biomart_web_file(self, tmp_path):
        """Create a BioMart export using the web interface's display labels."""
        content = """Chromosome/scaffold name\tGene start (bp)\tGene end (bp)\tGene name\tStrand
1\t1000000\t1020000\tBRCA1\t1
1\t1050000\t1080000\tTP53\t-1
"""
        filepath = tmp_path / "biomart_web.tsv"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def biomart_attr_file(self, tmp_path):
        """Create a BioMart export using the biomaRt attribute labels."""
        content = """chromosome_name\tstart_position\tend_position\texternal_gene_name\tstrand
1\t1000000\t1020000\tBRCA1\t1
1\t1050000\t1080000\tTP53\t-1
"""
        filepath = tmp_path / "biomart_attributes.tsv"
        filepath.write_text(content)
        return filepath

    def test_load_ensembl_genes_web_export_labels(self, biomart_web_file):
        """Test that BioMart display labels map to standard column names."""
        df = load_ensembl_genes(biomart_web_file)

        assert len(df) == 2
        assert df["chr"].iloc[0] == 1
        assert df["start"].iloc[0] == 1000000
        assert df["end"].iloc[0] == 1020000
        assert df["gene_name"].iloc[0] == "BRCA1"

    def test_load_ensembl_genes_web_export_strand_symbols(self, biomart_web_file):
        """Test that display-label strand integers become +/- symbols."""
        df = load_ensembl_genes(biomart_web_file)

        assert df["strand"].iloc[0] == "+"  # Not the raw 1
        assert df["strand"].iloc[1] == "-"  # Not the raw -1

    def test_load_ensembl_genes_attribute_labels(self, biomart_attr_file):
        """Test that biomaRt attribute labels map to standard column names."""
        df = load_ensembl_genes(biomart_attr_file)

        assert len(df) == 2
        assert df["chr"].iloc[0] == 1
        assert df["start"].iloc[0] == 1000000
        assert df["end"].iloc[0] == 1020000
        assert df["gene_name"].iloc[0] == "BRCA1"

    def test_load_ensembl_genes_attribute_strand_symbols(self, biomart_attr_file):
        """Test that attribute-label strand integers become +/- symbols."""
        df = load_ensembl_genes(biomart_attr_file)

        assert df["strand"].iloc[0] == "+"  # Not the raw 1
        assert df["strand"].iloc[1] == "-"  # Not the raw -1


# =============================================================================
# Test Validation Functions
# =============================================================================


class TestGWASValidation:
    """Tests for GWAS DataFrame validation."""

    def test_valid_gwas_df_passes(self):
        """Test that valid GWAS data passes validation."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000, 1002000],
                "p_wald": [0.01, 0.001, 1e-8],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )

        result = validate_gwas_dataframe(df)
        assert result is not None

    def test_missing_position_column_fails(self):
        """Test that missing position column raises error."""
        df = pd.DataFrame(
            {
                "p_wald": [0.01, 0.001],
                "rs": ["rs1", "rs2"],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing columns"):
            validate_gwas_dataframe(df)

    def test_missing_pvalue_column_fails(self):
        """Test that missing p-value column raises error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "rs": ["rs1", "rs2"],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing columns"):
            validate_gwas_dataframe(df)

    def test_negative_position_fails(self):
        """Test that negative positions raise error."""
        df = pd.DataFrame(
            {
                "ps": [-1000, 1001000],
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values <= 0"):
            validate_gwas_dataframe(df)

    def test_pvalue_out_of_range_fails(self):
        """Test that p-values outside (0, 1] raise error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": [0.01, 1.5],  # 1.5 is out of range
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values > 1"):
            validate_gwas_dataframe(df)

    def test_zero_pvalue_fails(self):
        """Test that p-value of 0 raises error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": [0.0, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values <= 0"):
            validate_gwas_dataframe(df)

    def test_nan_position_fails(self):
        """Test that NaN positions raise error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, None],
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"null values"):
            validate_gwas_dataframe(df)

    def test_non_numeric_pvalue_fails(self):
        """Test that non-numeric p-values raise clear validation error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": ["0.01", "significant"],  # Strings, not numbers
            }
        )

        with pytest.raises(LoaderValidationError, match="must be numeric"):
            validate_gwas_dataframe(df)

    def test_non_numeric_position_fails(self):
        """Test that non-numeric positions raise clear validation error."""
        df = pd.DataFrame(
            {
                "ps": ["chr1:1000", "chr1:2000"],  # Strings, not numbers
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match="must be numeric"):
            validate_gwas_dataframe(df)


class TestEQTLValidation:
    """Tests for eQTL DataFrame validation."""

    def test_valid_eqtl_df_passes(self):
        """Test that valid eQTL data passes validation."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
                "gene": ["BRCA1", "BRCA1"],
                "effect": [0.5, -0.3],
            }
        )

        result = validate_eqtl_dataframe(df)
        assert result is not None

    def test_missing_gene_column_fails(self):
        """Test that missing gene column raises error."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "p_value": [1e-6, 0.01],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing columns"):
            validate_eqtl_dataframe(df)


class TestFinemappingValidation:
    """Tests for fine-mapping DataFrame validation."""

    def test_valid_finemapping_df_passes(self):
        """Test that valid fine-mapping data passes validation."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000, 1002000],
                "pip": [0.85, 0.12, 0.03],
                "cs": [1, 1, 0],
            }
        )

        result = validate_finemapping_dataframe(df)
        assert result is not None

    def test_pip_out_of_range_fails(self):
        """Test that PIP outside [0, 1] raises error."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "pip": [0.85, 1.5],  # 1.5 is out of range
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values > 1"):
            validate_finemapping_dataframe(df)

    def test_negative_pip_fails(self):
        """Test that negative PIP raises error."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "pip": [-0.1, 0.5],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values < 0"):
            validate_finemapping_dataframe(df)


class TestGenesValidation:
    """Tests for genes DataFrame validation."""

    def test_valid_genes_df_passes(self):
        """Test that valid genes data passes validation."""
        df = pd.DataFrame(
            {
                "chr": ["1", "1", "1"],
                "start": [1000000, 1050000, 1100000],
                "end": [1020000, 1080000, 1150000],
                "gene_name": ["GENE1", "GENE2", "GENE3"],
            }
        )

        result = validate_genes_dataframe(df)
        assert result is not None

    def test_end_before_start_fails(self):
        """Test that end < start raises error."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1020000],  # Start after end
                "end": [1000000],
                "gene_name": ["GENE1"],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"rows have start > end"):
            validate_genes_dataframe(df)

    def test_negative_start_fails(self):
        """Test that negative start raises error."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [-1000],
                "end": [1000000],
                "gene_name": ["GENE1"],
            }
        )

        with pytest.raises(LoaderValidationError, match=r"values < 0"):
            validate_genes_dataframe(df)


# =============================================================================
# Test File Path Validation
# =============================================================================


class TestFileValidation:
    """Tests for file path validation."""

    def test_nonexistent_file_raises_error(self, tmp_path):
        """Test that non-existent file raises appropriate error."""
        fake_path = tmp_path / "nonexistent.assoc"

        with pytest.raises(Exception):  # FileNotFoundError or LoaderValidationError
            load_plink_assoc(fake_path)

    def test_corrupt_file_raises_error(self, tmp_path):
        """Test that corrupt/empty file raises an error during loading."""
        filepath = tmp_path / "corrupt.assoc"
        filepath.write_text("this is not valid tabular data\n!!!\n")

        with pytest.raises(Exception):
            load_plink_assoc(filepath)

    def test_empty_file_raises_error(self, tmp_path):
        """Test that empty file raises an error."""
        filepath = tmp_path / "empty.assoc"
        filepath.write_text("")

        with pytest.raises(Exception):
            load_plink_assoc(filepath)


# =============================================================================
# Integration Tests
# =============================================================================


class TestLoaderIntegration:
    """Integration tests for loader -> validation -> plotting flow."""

    def test_loaded_gwas_ready_for_plotting(self, plink_assoc_file):
        """Test that loaded GWAS data is ready for plotting."""
        df = load_plink_assoc(plink_assoc_file)

        # Should have required columns with correct types
        assert df["ps"].dtype in ["int64", "int32", "float64"]
        assert df["p_wald"].dtype == "float64"

        # Values should be in valid ranges
        assert (df["ps"] > 0).all()
        assert (df["p_wald"] > 0).all()
        assert (df["p_wald"] <= 1).all()

    def test_loaded_finemapping_ready_for_plotting(self, susie_file):
        """Test that loaded fine-mapping data is ready for plotting."""
        df = load_susie(susie_file)

        # Should have required columns
        assert "pos" in df.columns
        assert "pip" in df.columns

        # Values in valid ranges
        assert (df["pos"] > 0).all()
        assert (df["pip"] >= 0).all()
        assert (df["pip"] <= 1).all()

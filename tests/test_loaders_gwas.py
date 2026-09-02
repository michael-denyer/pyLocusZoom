"""Tests for the GWAS summary-statistic file loaders."""

import pytest

from pylocuszoom.loaders import (
    load_bolt_lmm,
    load_gemma,
    load_gwas,
    load_gwas_catalog,
    load_plink_assoc,
    load_regenie,
    load_saige,
)


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
        assert df["pos"].iloc[0] == 1000000
        assert df["rs"].iloc[0] == "rs123"
        assert df["p_value"].iloc[0] == 0.01

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
        assert df["pos"].iloc[0] == 1000000  # POS, the second pos_col candidate
        assert df["rs"].iloc[0] == "rs123"  # ID, the second rs_col candidate
        assert df["chr"].iloc[0] == 1  # CHROM, the third chr candidate
        assert df["p_value"].iloc[0] == 0.01

    def test_load_plink_assoc_basic(self, plink_assoc_file):
        """Test basic PLINK file loading."""
        df = load_plink_assoc(plink_assoc_file)

        assert "pos" in df.columns
        assert "p_value" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_plink_assoc_custom_columns(self, plink_assoc_file):
        """Naming the output columns still works, and warns that it will not."""
        with pytest.warns(DeprecationWarning, match="5.0.0"):
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

        assert df["pos"].iloc[0] == 1000000
        assert df["p_value"].iloc[0] == 0.01
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

        assert df["p_value"].iloc[0] == 0.004  # Taken as-is, not 10**(-0.004)
        assert df["p_value"].iloc[2] == pytest.approx(2e-9, rel=0.01)

    def test_load_regenie_basic(self, regenie_file):
        """Test basic REGENIE file loading."""
        df = load_regenie(regenie_file)

        assert "pos" in df.columns
        assert "p_value" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_regenie_log10p_conversion(self, regenie_file):
        """Test that LOG10P is converted to p-value."""
        df = load_regenie(regenie_file)

        # LOG10P=2.0 -> p=0.01, LOG10P=8.0 -> p=1e-8
        assert df["p_value"].iloc[0] == pytest.approx(0.01, rel=0.01)
        assert df["p_value"].iloc[2] == pytest.approx(1e-8, rel=0.01)


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

        assert "p_value" in df.columns
        # Should use SPA-adjusted p-values (p.value.NA column)
        assert df["p_value"].iloc[0] == 0.005  # Not 0.01
        assert df["p_value"].iloc[2] == pytest.approx(1e-9, rel=0.01)

    def test_load_saige_fallback_to_pvalue(self, saige_file_no_spa):
        """Test that SAIGE loader falls back to p.value when p.value.NA missing."""
        df = load_saige(saige_file_no_spa)

        assert "p_value" in df.columns
        assert df["p_value"].iloc[0] == 0.01

    def test_load_saige_basic(self, saige_file):
        """Test basic SAIGE file loading."""
        df = load_saige(saige_file)

        assert "pos" in df.columns
        assert "p_value" in df.columns
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

        assert "pos" in df.columns
        assert "p_value" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_bolt_prefers_full_model(self, bolt_file):
        """Test that BOLT-LMM loader prefers P_BOLT_LMM over P_BOLT_LMM_INF."""
        df = load_bolt_lmm(bolt_file)

        # Should use P_BOLT_LMM (full model), not P_BOLT_LMM_INF
        assert df["p_value"].iloc[0] == 0.005  # Not 0.01
        assert df["p_value"].iloc[2] == pytest.approx(1e-9, rel=0.01)

    def test_load_bolt_fallback_to_inf(self, bolt_file_inf_only):
        """Test that BOLT-LMM loader falls back to P_BOLT_LMM_INF."""
        df = load_bolt_lmm(bolt_file_inf_only)

        assert "p_value" in df.columns
        assert df["p_value"].iloc[0] == 0.01


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

        assert df["p_value"].iloc[0] == 0.001  # Not 0.02 (p_lrt), not 0.3 (p_score)

    def test_load_gemma_precedence_with_custom_p_col(self, gemma_all_p_file):
        """Test GEMMA p-value precedence through a custom output column.

        With the default p_col a reordered candidate tuple would rename p_lrt
        onto the existing p_wald column, so the assertion would hit duplicate
        labels rather than a clean value mismatch.
        """
        with pytest.warns(DeprecationWarning):
            df = load_gemma(gemma_all_p_file, p_col="pval")

        assert df["pval"].iloc[0] == 0.001  # Not 0.02 (p_lrt), not 0.3 (p_score)

    def test_load_gemma_basic(self, gemma_file):
        """Test basic GEMMA file loading."""
        df = load_gemma(gemma_file)

        assert "pos" in df.columns
        assert "p_value" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_gemma_values(self, gemma_file):
        """Test that values are loaded correctly."""
        df = load_gemma(gemma_file)

        assert df["pos"].iloc[0] == 1000000
        assert df["p_value"].iloc[0] == 0.01

    def test_load_gemma_fallback_to_lrt(self, gemma_lrt_file):
        """Test GEMMA loader falls back to p_lrt when p_wald missing."""
        df = load_gemma(gemma_lrt_file)

        assert "p_value" in df.columns
        assert df["p_value"].iloc[0] == 0.02


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

        assert "pos" in df.columns
        assert "p_value" in df.columns
        assert "rs" in df.columns
        assert len(df) == 3

    def test_load_gwas_catalog_values(self, catalog_file):
        """Test that values are mapped correctly."""
        df = load_gwas_catalog(catalog_file)

        assert df["pos"].iloc[0] == 1000000
        assert df["p_value"].iloc[0] == 0.01
        assert df["rs"].iloc[0] == "rs123"

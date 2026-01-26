"""Tests for file format loaders and pydantic validation."""

import pandas as pd
import pytest

from pylocuszoom.loaders import (
    load_bed,
    load_gwas,
    load_plink_assoc,
    load_regenie,
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


# =============================================================================
# Test Fine-mapping Loaders
# =============================================================================


class TestSuSiELoader:
    """Tests for SuSiE file loader."""

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

        with pytest.raises(LoaderValidationError, match="Missing required column"):
            validate_gwas_dataframe(df)

    def test_missing_pvalue_column_fails(self):
        """Test that missing p-value column raises error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "rs": ["rs1", "rs2"],
            }
        )

        with pytest.raises(LoaderValidationError, match="Missing required column"):
            validate_gwas_dataframe(df)

    def test_negative_position_fails(self):
        """Test that negative positions raise error."""
        df = pd.DataFrame(
            {
                "ps": [-1000, 1001000],
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match="non-positive"):
            validate_gwas_dataframe(df)

    def test_pvalue_out_of_range_fails(self):
        """Test that p-values outside (0, 1] raise error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": [0.01, 1.5],  # 1.5 is out of range
            }
        )

        with pytest.raises(LoaderValidationError, match="outside range"):
            validate_gwas_dataframe(df)

    def test_zero_pvalue_fails(self):
        """Test that p-value of 0 raises error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, 1001000],
                "p_wald": [0.0, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match="outside range"):
            validate_gwas_dataframe(df)

    def test_nan_position_fails(self):
        """Test that NaN positions raise error."""
        df = pd.DataFrame(
            {
                "ps": [1000000, None],
                "p_wald": [0.01, 0.001],
            }
        )

        with pytest.raises(LoaderValidationError, match="missing values"):
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

        with pytest.raises(LoaderValidationError, match="Missing required column"):
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

        with pytest.raises(LoaderValidationError, match="outside range"):
            validate_finemapping_dataframe(df)

    def test_negative_pip_fails(self):
        """Test that negative PIP raises error."""
        df = pd.DataFrame(
            {
                "pos": [1000000, 1001000],
                "pip": [-0.1, 0.5],
            }
        )

        with pytest.raises(LoaderValidationError, match="outside range"):
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

        with pytest.raises(LoaderValidationError, match="end < start"):
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

        with pytest.raises(LoaderValidationError, match="negative"):
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

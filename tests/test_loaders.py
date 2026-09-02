"""Tests for loader dispatch, format detection and file-path validation.

Per-family loader tests live in ``test_loaders_gwas.py``,
``test_loaders_eqtl.py``, ``test_loaders_finemapping.py`` and
``test_loaders_annotation.py``.
"""

import io

import pytest

from pylocuszoom.loaders import (
    load_gwas,
    load_plink_assoc,
    load_susie,
)
from pylocuszoom.loaders.gwas import _detect_format


class TestAutoFormatDetection:
    """Tests for automatic format detection."""

    def test_load_gwas_detects_plink(self, plink_assoc_file):
        """Test that load_gwas auto-detects PLINK format."""
        df = load_gwas(plink_assoc_file)
        assert "pos" in df.columns
        assert len(df) == 3

    def test_load_gwas_detects_regenie(self, regenie_file):
        """Test that load_gwas auto-detects REGENIE format."""
        df = load_gwas(regenie_file)
        assert "pos" in df.columns
        assert len(df) == 3

    def test_load_gwas_detects_bolt_from_stats_ext(self, tmp_path):
        """Test that .stats extension is detected as BOLT-LMM format."""
        content = """SNP\tCHR\tBP\tALLELE1\tALLELE0\tA1FREQ\tBETA\tSE\tP_BOLT_LMM_INF
rs123\t1\t1000000\tA\tG\t0.3\t0.5\t0.2\t0.01
"""
        filepath = tmp_path / "result.stats"
        filepath.write_text(content)
        df = load_gwas(filepath)
        assert "pos" in df.columns
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

        assert "pos" in df.columns
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
        assert "pos" in df.columns

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
        assert "pos" in df.columns

    def test_load_gwas_detects_qassoc(self, tmp_path):
        """Test .qassoc extension is detected as plink."""
        content = """CHR SNP BP A1 TEST NMISS BETA STAT P
1 rs123 1000000 A ADD 1000 0.5 2.5 0.01
"""
        filepath = tmp_path / "result.qassoc"
        filepath.write_text(content)
        df = load_gwas(filepath)
        assert "pos" in df.columns

    def test_load_gwas_detects_gemma_from_assoc_txt(self, tmp_path):
        """A GEMMA .assoc.txt file must not be claimed by the shorter .assoc hint."""
        content = """chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tlogl_H1\tl_remle\tp_wald
1\trs123\t1000000\t0\tA\tG\t0.3\t0.5\t0.2\t100\t0.5\t0.01
"""
        filepath = tmp_path / "output.assoc.txt"
        filepath.write_text(content)

        assert _detect_format(filepath) == "gemma"
        assert _detect_format(tmp_path / "result.assoc") == "plink"
        df = load_gwas(filepath)
        assert "pos" in df.columns
        assert len(df) == 1


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


class TestLoaderIntegration:
    """Integration tests for loader -> validation -> plotting flow."""

    def test_loaded_gwas_ready_for_plotting(self, plink_assoc_file):
        """Test that loaded GWAS data is ready for plotting."""
        df = load_plink_assoc(plink_assoc_file)

        # Should have required columns with correct types
        assert df["pos"].dtype in ["int64", "int32", "float64"]
        assert df["p_value"].dtype == "float64"

        # Values should be in valid ranges
        assert (df["pos"] > 0).all()
        assert (df["p_value"] > 0).all()
        assert (df["p_value"] <= 1).all()

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

"""Tests for LD calculation module."""

import subprocess
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from pylocuszoom.exceptions import EmptyLDOutputError, PlinkError, ValidationError
from pylocuszoom.ld import (
    _add_species_flags,
    build_ld_command,
    build_pairwise_ld_command,
    calculate_ld,
    calculate_pairwise_ld,
    find_plink,
    parse_ld_output,
    parse_pairwise_ld_output,
    validate_plink_files,
)


class TestValidatePlinkFiles:
    """Tests for validate_plink_files function."""

    def test_valid_plink_files(self, tmp_path):
        """Valid PLINK fileset passes."""
        # Create all required files
        (tmp_path / "test.bed").touch()
        (tmp_path / "test.bim").touch()
        (tmp_path / "test.fam").touch()

        result = validate_plink_files(tmp_path / "test")
        assert result == tmp_path / "test"

    def test_missing_bed_raises(self, tmp_path):
        """Missing .bed file raises error."""
        (tmp_path / "test.bim").touch()
        (tmp_path / "test.fam").touch()

        with pytest.raises(ValidationError, match=".bed"):
            validate_plink_files(tmp_path / "test")

    def test_missing_multiple_files(self, tmp_path):
        """Missing multiple files lists all in error."""
        (tmp_path / "test.bed").touch()
        # Missing .bim and .fam

        with pytest.raises(ValidationError) as exc_info:
            validate_plink_files(tmp_path / "test")

        assert ".bim" in str(exc_info.value)
        assert ".fam" in str(exc_info.value)

    def test_prefix_with_dots_preserved(self, tmp_path):
        """Prefixes containing dots (e.g. 'ukbb.v3') must not be truncated.

        Regression: an earlier implementation used Path.with_suffix(),
        which would rewrite 'ukbb.v3' -> 'ukbb.bed', checking the wrong
        file on disk. Real files with the full prefix would appear
        missing, and the existence check would pass only against files
        that don't exist.
        """
        prefix = tmp_path / "ukbb.v3"
        (tmp_path / "ukbb.v3.bed").touch()
        (tmp_path / "ukbb.v3.bim").touch()
        (tmp_path / "ukbb.v3.fam").touch()

        result = validate_plink_files(prefix)
        assert result == prefix

    def test_prefix_with_dots_missing_raises(self, tmp_path):
        """Dot-containing prefix with missing files raises (not silently passes)."""
        prefix = tmp_path / "ukbb.v3"
        (tmp_path / "ukbb.bed").touch()
        (tmp_path / "ukbb.bim").touch()
        (tmp_path / "ukbb.fam").touch()

        with pytest.raises(ValidationError):
            validate_plink_files(prefix)


class TestFindPlink:
    """Tests for find_plink function."""

    def test_returns_plink_path_when_found(self):
        """Should return path when PLINK is on PATH."""
        with patch("shutil.which") as mock_which:
            mock_which.return_value = "/usr/bin/plink1.9"
            result = find_plink()
            assert result == "/usr/bin/plink1.9"

    def test_tries_plink19_first(self):
        """Should try plink1.9 before plink."""
        with patch("shutil.which") as mock_which:
            mock_which.side_effect = lambda x: (
                "/usr/bin/plink1.9" if x == "plink1.9" else None
            )
            result = find_plink()
            assert result == "/usr/bin/plink1.9"

    def test_falls_back_to_plink(self):
        """Should fall back to plink if plink1.9 not found."""
        with patch("shutil.which") as mock_which:
            mock_which.side_effect = lambda x: (
                "/usr/bin/plink" if x == "plink" else None
            )
            result = find_plink()
            assert result == "/usr/bin/plink"

    def test_returns_none_when_not_found(self):
        """Should return None when PLINK not on PATH."""
        with patch("shutil.which", return_value=None):
            result = find_plink()
            assert result is None


class TestBuildLdCommand:
    """Tests for build_ld_command function."""

    def test_includes_required_flags(self):
        """Command should include all required PLINK flags."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
        )

        assert "/usr/bin/plink1.9" in cmd
        assert "--bfile" in cmd
        assert "/path/to/data" in cmd
        assert "--r2" in cmd
        assert "--ld-snp" in cmd
        assert "rs12345" in cmd
        assert "--out" in cmd
        assert "/path/to/output" in cmd

    def test_window_kb_parameter(self):
        """Command should include specified window size."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            window_kb=1000,
        )
        assert "--ld-window-kb" in cmd
        idx = cmd.index("--ld-window-kb")
        assert cmd[idx + 1] == "1000"

    def test_removes_default_snp_limit(self):
        """Command should set --ld-window 99999 to remove 10 SNP default."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
        )
        assert "--ld-window" in cmd
        idx = cmd.index("--ld-window")
        assert cmd[idx + 1] == "99999"

    def test_includes_threads(self):
        """Command should include thread count."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            threads=4,
        )
        assert "--threads" in cmd
        idx = cmd.index("--threads")
        assert cmd[idx + 1] == "4"


class TestParseLdOutput:
    """Tests for parse_ld_output function."""

    def test_parses_plink_whitespace_separated_output(self, tmp_path):
        """Should parse PLINK's whitespace-separated .ld file."""
        ld_content = """CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       1000    rs12345 1       1500    rs11111 0.95
1       1000    rs12345 1       2000    rs22222 0.75
1       1000    rs12345 1       2500    rs33333 0.45"""

        ld_file = tmp_path / "test.ld"
        ld_file.write_text(ld_content)

        result = parse_ld_output(str(ld_file), "rs12345")

        assert len(result) == 4  # 3 SNPs + lead SNP
        assert "SNP" in result.columns
        assert "R2" in result.columns

        # Check parsed values
        snps = result["SNP"].tolist()
        assert "rs11111" in snps
        assert "rs22222" in snps
        assert "rs33333" in snps
        assert "rs12345" in snps  # Lead SNP added

    def test_adds_lead_snp_with_r2_one(self, tmp_path):
        """Should add lead SNP with R2=1.0."""
        ld_content = """CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       1000    rs12345 1       1500    rs11111 0.95"""

        ld_file = tmp_path / "test.ld"
        ld_file.write_text(ld_content)

        result = parse_ld_output(str(ld_file), "rs12345")

        lead_row = result[result["SNP"] == "rs12345"]
        assert len(lead_row) == 1
        assert lead_row["R2"].iloc[0] == 1.0

    def test_deduplicates_lead_self_pair(self, tmp_path):
        """PLINK's --ld-snp output includes the lead paired with itself (R2=1.0).

        Regression: parse_ld_output appended an explicit lead row on top of that
        self-pair, so the lead SNP appeared twice. The duplicate key later broke
        the ``validate="many_to_one"`` LD merge in ``LocusZoomPlotter.plot()``
        with a ``MergeError``, silently disabling LD colouring for every real
        PLINK run (plink1.9 emits the self-pair; the older fixtures did not).
        """
        ld_content = """CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       1000    rs12345 1       1000    rs12345 1.0
1       1000    rs12345 1       1500    rs11111 0.95
1       1000    rs12345 1       2000    rs22222 0.75"""

        ld_file = tmp_path / "selfpair.ld"
        ld_file.write_text(ld_content)

        result = parse_ld_output(str(ld_file), "rs12345")

        assert result["SNP"].is_unique  # no duplicate merge key for the lead
        lead_row = result[result["SNP"] == "rs12345"]
        assert len(lead_row) == 1
        assert lead_row["R2"].iloc[0] == 1.0

    def test_raises_for_missing_file(self, tmp_path):
        """Should raise PlinkError for missing output file."""
        with pytest.raises(PlinkError, match="output file not found"):
            parse_ld_output(str(tmp_path / "nonexistent.ld"), "rs12345")

    def test_raises_for_empty_ld_output(self, tmp_path):
        """Empty .ld file (header only) must raise, not return empty DataFrame.

        Regression: PLINK exits 0 but writes no LD pairs when the lead SNP
        is monomorphic, filtered by --maf, or absent from the reference.
        Silently returning an empty DataFrame left callers with an
        uncoloured plot and no diagnostic.
        """
        ld_content = "CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2\n"
        ld_file = tmp_path / "empty.ld"
        ld_file.write_text(ld_content)

        with pytest.raises(EmptyLDOutputError, match="empty LD output"):
            parse_ld_output(str(ld_file), "rs12345")

    def test_parses_r2_boundary_values(self, tmp_path):
        """Should correctly parse R2 boundary values."""
        ld_content = """CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       1000    rs12345 1       1500    rs11111 1.0
1       1000    rs12345 1       2000    rs22222 0.0
1       1000    rs12345 1       2500    rs33333 0.5"""

        ld_file = tmp_path / "test.ld"
        ld_file.write_text(ld_content)

        result = parse_ld_output(str(ld_file), "rs12345")

        r2_values = result[result["SNP"] != "rs12345"]["R2"].tolist()
        assert 1.0 in r2_values
        assert 0.0 in r2_values
        assert 0.5 in r2_values


class TestCalculateLd:
    """Tests for calculate_ld function."""

    @pytest.fixture
    def mock_plink_files(self, tmp_path):
        """Create mock PLINK files for testing."""
        bfile = tmp_path / "test_geno"
        # Create empty placeholder files
        (bfile.parent / f"{bfile.name}.bed").touch()
        (bfile.parent / f"{bfile.name}.bim").touch()
        (bfile.parent / f"{bfile.name}.fam").touch()
        return str(bfile)

    def test_parses_the_ld_body_plink_wrote(self, tmp_path, fake_plink):
        """A successful run returns the R² frame parsed out of PLINK's own file."""
        bfile, plink_writes = fake_plink
        body = """CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       1000    rs12345 1       1000    rs12345 1.0
1       1000    rs12345 1       1500    rs11111 0.95
1       1000    rs12345 1       2000    rs22222 0.40"""

        with plink_writes(body):
            result = calculate_ld(
                bfile_path=bfile,
                lead_snp="rs12345",
                plink_path="/usr/bin/plink1.9",
                working_dir=str(tmp_path),
            )

        assert list(result.columns) == ["SNP", "R2"]
        assert dict(zip(result["SNP"], result["R2"], strict=True)) == {
            "rs12345": 1.0,
            "rs11111": 0.95,
            "rs22222": 0.40,
        }

    def test_a_colon_separated_snp_id_survives_the_round_trip(
        self, tmp_path, fake_plink
    ):
        """A VCF-style lead id names a file PLINK can write and the parser can read.

        The id goes into the output path, so a colon has to be sanitised
        there while the id PLINK is asked for, and the one the parser looks
        up in the result, stay verbatim.
        """
        bfile, plink_writes = fake_plink
        body = """CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       12345   1:12345:A:G     1       1500    rs11111 0.80"""

        with plink_writes(body):
            result = calculate_ld(
                bfile_path=bfile,
                lead_snp="1:12345:A:G",
                plink_path="/usr/bin/plink1.9",
                working_dir=str(tmp_path),
            )

        lead = result[result["SNP"] == "1:12345:A:G"]
        assert len(lead) == 1
        assert lead["R2"].iloc[0] == 1.0

    def test_an_empty_ld_body_raises_rather_than_colouring_nothing(
        self, tmp_path, fake_plink
    ):
        """PLINK exits 0 with a header-only file when the lead is monomorphic."""
        bfile, plink_writes = fake_plink

        with plink_writes("CHR_A   BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2\n"):
            with pytest.raises(EmptyLDOutputError, match="empty LD output"):
                calculate_ld(
                    bfile_path=bfile,
                    lead_snp="rs12345",
                    plink_path="/usr/bin/plink1.9",
                    working_dir=str(tmp_path),
                )

    def test_raises_when_plink_not_found(self, mock_plink_files):
        """Should raise FileNotFoundError when PLINK not found."""
        with patch("pylocuszoom.ld.find_plink", return_value=None):
            with pytest.raises(FileNotFoundError, match="PLINK not found"):
                calculate_ld(
                    bfile_path=mock_plink_files,
                    lead_snp="rs12345",
                )

    def test_raises_plink_error_on_plink_failure(self, tmp_path, fake_plink):
        """Should raise PlinkError when PLINK returns non-zero exit code."""
        bfile, plink_writes = fake_plink

        with plink_writes(None, returncode=1, stderr="Error: invalid SNP"):
            with pytest.raises(PlinkError, match="exit code 1"):
                calculate_ld(
                    bfile_path=bfile,
                    lead_snp="rs12345",
                    plink_path="/usr/bin/plink1.9",
                    working_dir=str(tmp_path),
                )

    def test_plink_error_includes_stderr_in_message(self, tmp_path, fake_plink):
        """PlinkError message should include PLINK's stderr output."""
        bfile, plink_writes = fake_plink

        with plink_writes(None, returncode=1, stderr="Error: variant not found"):
            with pytest.raises(PlinkError, match="variant not found"):
                calculate_ld(
                    bfile_path=bfile,
                    lead_snp="rs12345",
                    plink_path="/usr/bin/plink1.9",
                    working_dir=str(tmp_path),
                )

    def test_raises_plink_error_on_timeout(self, tmp_path, mock_plink_files):
        """Should raise PlinkError when PLINK times out."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.side_effect = subprocess.TimeoutExpired(
                    cmd="plink", timeout=300
                )

                with pytest.raises(PlinkError, match="timed out"):
                    calculate_ld(
                        bfile_path=mock_plink_files,
                        lead_snp="rs12345",
                        working_dir=str(tmp_path),
                    )

    def test_cleans_up_temp_directory(self, tmp_path, monkeypatch, mock_plink_files):
        """A failed run leaves nothing behind in the directory it created."""
        temp_base = tmp_path / "tmpbase"
        temp_base.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(temp_base))

        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                with pytest.raises(PlinkError):
                    calculate_ld(
                        bfile_path=mock_plink_files,
                        lead_snp="rs12345",
                        working_dir=None,
                    )

        assert list(temp_base.iterdir()) == []

    def test_raises_validation_error_for_missing_plink_files(self, tmp_path):
        """Bug: calculate_ld() raises ValidationError for missing PLINK files.

        The docstring only documents FileNotFoundError, but validate_plink_files()
        raises ValidationError when .bed/.bim/.fam files are missing.
        This test documents the actual behavior.
        """
        from pylocuszoom.utils import ValidationError

        # Non-existent PLINK files
        nonexistent_bfile = str(tmp_path / "nonexistent")

        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            # Should raise ValidationError (not FileNotFoundError as docstring says)
            with pytest.raises(ValidationError, match="PLINK files missing"):
                calculate_ld(
                    bfile_path=nonexistent_bfile,
                    lead_snp="rs12345",
                )


class TestBuildPairwiseLdCommand:
    """Tests for build_pairwise_ld_command function."""

    def test_includes_r2_square_flag(self):
        """Command should include --r2 square for pairwise matrix."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
        )

        assert "--r2" in cmd
        assert "square" in cmd

    def test_includes_write_snplist_flag(self):
        """Command should include --write-snplist to track SNP order."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
        )

        assert "--write-snplist" in cmd

    def test_includes_extract_with_snp_list_file(self):
        """Command should include --extract when snp_list_file provided."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            snp_list_file="/path/to/snps.txt",
        )

        assert "--extract" in cmd
        assert "/path/to/snps.txt" in cmd

    def test_includes_region_flags_when_provided(self):
        """Command should include --chr, --from-bp, --to-bp for region mode."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            chrom=1,
            start=1000000,
            end=2000000,
        )

        assert "--chr" in cmd
        assert "1" in cmd
        assert "--from-bp" in cmd
        assert "1000000" in cmd
        assert "--to-bp" in cmd
        assert "2000000" in cmd

    def test_uses_dprime_metric(self):
        """Command should use --r dprime for D' metric."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            metric="dprime",
        )

        # For D', PLINK uses --r (not --r2) with dprime modifier
        assert "--r" in cmd
        assert "dprime" in cmd
        assert "square" in cmd

    @pytest.mark.parametrize("metric", ["Dprime", "r", "R2", ""])
    def test_unknown_metric_is_rejected(self, metric):
        """Any spelling other than the two accepted metrics is an error, not r2."""
        with pytest.raises(ValidationError, match="dprime"):
            build_pairwise_ld_command(
                plink_path="/usr/bin/plink1.9",
                bfile_path="/path/to/data",
                output_path="/path/to/output",
                metric=metric,
            )

    def test_calculate_pairwise_ld_rejects_metric_before_running_plink(self, tmp_path):
        """A bad metric fails at the boundary, before PLINK is even looked up."""
        for ext in (".bed", ".bim", ".fam"):
            (tmp_path / f"data{ext}").touch()
        with patch("pylocuszoom.ld.subprocess.run") as mock_run:
            with pytest.raises(ValidationError, match="dprime"):
                calculate_pairwise_ld(
                    bfile_path=str(tmp_path / "data"),
                    plink_path="/usr/bin/plink1.9",
                    metric="Dprime",
                )
        mock_run.assert_not_called()

    def test_default_metric_is_r2(self):
        """Command should use --r2 by default."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
        )

        assert "--r2" in cmd


class TestParsePairwiseLdOutput:
    """Tests for parse_pairwise_ld_output function."""

    def test_parses_square_matrix_format(self, tmp_path):
        """Should parse PLINK's square matrix .ld file (no header, whitespace-separated)."""
        # PLINK --r2 square outputs N x N matrix without headers
        ld_content = """1.0\t0.85\t0.45
0.85\t1.0\t0.60
0.45\t0.60\t1.0"""
        snplist_content = """rs1
rs2
rs3"""

        ld_file = tmp_path / "test.ld"
        snplist_file = tmp_path / "test.snplist"
        ld_file.write_text(ld_content)
        snplist_file.write_text(snplist_content)

        matrix, snp_ids = parse_pairwise_ld_output(str(ld_file), str(snplist_file))

        assert matrix.shape == (3, 3)
        assert snp_ids == ["rs1", "rs2", "rs3"]

    def test_matrix_has_snp_ids_as_index_and_columns(self, tmp_path):
        """Matrix should have SNP IDs as both row index and column names."""
        ld_content = """1.0\t0.85
0.85\t1.0"""
        snplist_content = """rs1
rs2"""

        ld_file = tmp_path / "test.ld"
        snplist_file = tmp_path / "test.snplist"
        ld_file.write_text(ld_content)
        snplist_file.write_text(snplist_content)

        matrix, snp_ids = parse_pairwise_ld_output(str(ld_file), str(snplist_file))

        assert list(matrix.index) == ["rs1", "rs2"]
        assert list(matrix.columns) == ["rs1", "rs2"]

    def test_matrix_diagonal_is_one(self, tmp_path):
        """Diagonal elements (self-LD) should all be 1.0."""
        ld_content = """1.0\t0.85\t0.45
0.85\t1.0\t0.60
0.45\t0.60\t1.0"""
        snplist_content = """rs1
rs2
rs3"""

        ld_file = tmp_path / "test.ld"
        snplist_file = tmp_path / "test.snplist"
        ld_file.write_text(ld_content)
        snplist_file.write_text(snplist_content)

        matrix, _ = parse_pairwise_ld_output(str(ld_file), str(snplist_file))

        for snp in ["rs1", "rs2", "rs3"]:
            assert matrix.loc[snp, snp] == 1.0

    def test_matrix_is_symmetric(self, tmp_path):
        """Matrix should be symmetric (LD(A,B) == LD(B,A))."""
        ld_content = """1.0\t0.85\t0.45
0.85\t1.0\t0.60
0.45\t0.60\t1.0"""
        snplist_content = """rs1
rs2
rs3"""

        ld_file = tmp_path / "test.ld"
        snplist_file = tmp_path / "test.snplist"
        ld_file.write_text(ld_content)
        snplist_file.write_text(snplist_content)

        matrix, _ = parse_pairwise_ld_output(str(ld_file), str(snplist_file))

        # Check symmetry
        assert matrix.loc["rs1", "rs2"] == matrix.loc["rs2", "rs1"]
        assert matrix.loc["rs1", "rs3"] == matrix.loc["rs3", "rs1"]
        assert matrix.loc["rs2", "rs3"] == matrix.loc["rs3", "rs2"]

    def test_handles_nan_values(self, tmp_path):
        """Should handle nan values for SNP pairs without LD data."""
        ld_content = """1.0\tnan\t0.45
nan\t1.0\t0.60
0.45\t0.60\t1.0"""
        snplist_content = """rs1
rs2
rs3"""

        ld_file = tmp_path / "test.ld"
        snplist_file = tmp_path / "test.snplist"
        ld_file.write_text(ld_content)
        snplist_file.write_text(snplist_content)

        matrix, _ = parse_pairwise_ld_output(str(ld_file), str(snplist_file))

        import math

        assert math.isnan(matrix.loc["rs1", "rs2"])
        assert math.isnan(matrix.loc["rs2", "rs1"])

    def test_raises_for_missing_ld_file(self, tmp_path):
        """Should raise PlinkError for missing .ld file."""
        snplist_content = """rs1
rs2"""
        snplist_file = tmp_path / "test.snplist"
        snplist_file.write_text(snplist_content)

        with pytest.raises(PlinkError, match="output files missing"):
            parse_pairwise_ld_output(
                str(tmp_path / "nonexistent.ld"), str(snplist_file)
            )

    def test_raises_for_missing_snplist_file(self, tmp_path):
        """Should raise PlinkError for missing .snplist file."""
        ld_content = """1.0\t0.85
0.85\t1.0"""
        ld_file = tmp_path / "test.ld"
        ld_file.write_text(ld_content)

        with pytest.raises(PlinkError, match="output files missing"):
            parse_pairwise_ld_output(
                str(ld_file), str(tmp_path / "nonexistent.snplist")
            )

    def test_parses_space_separated_values(self, tmp_path):
        """Should parse space-separated values (not just tab-separated)."""
        ld_content = """1.0 0.85 0.45
0.85 1.0 0.60
0.45 0.60 1.0"""
        snplist_content = """rs1
rs2
rs3"""

        ld_file = tmp_path / "test.ld"
        snplist_file = tmp_path / "test.snplist"
        ld_file.write_text(ld_content)
        snplist_file.write_text(snplist_content)

        matrix, snp_ids = parse_pairwise_ld_output(str(ld_file), str(snplist_file))

        assert matrix.shape == (3, 3)
        assert matrix.loc["rs1", "rs2"] == 0.85


class TestCalculatePairwiseLd:
    """Tests for calculate_pairwise_ld function."""

    @pytest.fixture
    def mock_plink_files(self, tmp_path):
        """Create mock PLINK files for testing."""
        bfile = tmp_path / "test_geno"
        (bfile.parent / f"{bfile.name}.bed").touch()
        (bfile.parent / f"{bfile.name}.bim").touch()
        (bfile.parent / f"{bfile.name}.fam").touch()
        return str(bfile)

    def test_raises_when_plink_not_found(self, mock_plink_files):
        """Should raise FileNotFoundError when PLINK not found."""
        with patch("pylocuszoom.ld.find_plink", return_value=None):
            with pytest.raises(FileNotFoundError, match="PLINK not found"):
                calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2"],
                )

    def test_raises_plink_error_on_plink_failure(self, tmp_path, mock_plink_files):
        """Should raise PlinkError when PLINK returns non-zero exit code."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                with pytest.raises(PlinkError, match="exit code 1"):
                    calculate_pairwise_ld(
                        bfile_path=mock_plink_files,
                        snp_list=["rs1", "rs2"],
                        working_dir=str(tmp_path),
                    )

    def test_writes_snp_list_file(self, tmp_path, mock_plink_files):
        """Should write SNP list to file when snp_list provided."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                with pytest.raises(PlinkError):
                    calculate_pairwise_ld(
                        bfile_path=mock_plink_files,
                        snp_list=["rs1", "rs2", "rs3"],
                        working_dir=str(tmp_path),
                    )

                # Check the SNP list file was written (before PLINK ran)
                snp_list_file = tmp_path / "snp_list.txt"
                assert snp_list_file.exists()
                content = snp_list_file.read_text()
                assert "rs1" in content
                assert "rs2" in content
                assert "rs3" in content

    def test_region_mode_returns_every_snp_plink_kept(self, tmp_path, mock_plink_files):
        """Given a region and no SNP list, the result is whatever PLINK retained.

        Region mode writes no snp_list.txt, so there is nothing to check the
        result against and no SNP can go missing.
        """
        (tmp_path / "pairwise_ld.ld").write_text("1.0\t0.70\n0.70\t1.0")
        (tmp_path / "pairwise_ld.snplist").write_text("rs1\nrs2")

        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stderr="")

                matrix, snp_ids = calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    chrom=1,
                    start=1000000,
                    end=2000000,
                    working_dir=str(tmp_path),
                )

        assert snp_ids == ["rs1", "rs2"]
        assert matrix.loc["rs1", "rs2"] == 0.70
        assert not (tmp_path / "snp_list.txt").exists()

    def test_raises_validation_error_for_missing_snps(self, tmp_path, mock_plink_files):
        """Should raise ValidationError when requested SNPs not in reference."""
        from pylocuszoom.utils import ValidationError

        # Create output files with only rs1 and rs2 (missing rs3)
        ld_file = tmp_path / "pairwise_ld.ld"
        snplist_file = tmp_path / "pairwise_ld.snplist"
        ld_file.write_text("1.0\t0.85\n0.85\t1.0")
        snplist_file.write_text("rs1\nrs2")

        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stderr="")

                with pytest.raises(
                    ValidationError, match="SNPs not found in reference panel"
                ):
                    calculate_pairwise_ld(
                        bfile_path=mock_plink_files,
                        snp_list=["rs1", "rs2", "rs3"],  # rs3 not in output
                        working_dir=str(tmp_path),
                    )

    def test_cleans_up_temp_directory(self, tmp_path, monkeypatch, mock_plink_files):
        """A failed run leaves nothing behind in the directory it created."""
        temp_base = tmp_path / "tmpbase"
        temp_base.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(temp_base))

        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                with pytest.raises(PlinkError):
                    calculate_pairwise_ld(
                        bfile_path=mock_plink_files,
                        snp_list=["rs1", "rs2"],
                        working_dir=None,
                    )

        assert list(temp_base.iterdir()) == []

    def test_returns_matrix_and_snp_ids_on_success(self, tmp_path, mock_plink_files):
        """Should return (matrix, snp_ids) tuple on successful computation."""
        # Create output files
        ld_file = tmp_path / "pairwise_ld.ld"
        snplist_file = tmp_path / "pairwise_ld.snplist"
        ld_file.write_text("1.0\t0.85\t0.45\n0.85\t1.0\t0.60\n0.45\t0.60\t1.0")
        snplist_file.write_text("rs1\nrs2\nrs3")

        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0, stderr="")

                matrix, snp_ids = calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2", "rs3"],
                    working_dir=str(tmp_path),
                )

                assert matrix.shape == (3, 3)
                assert snp_ids == ["rs1", "rs2", "rs3"]
                assert matrix.loc["rs1", "rs2"] == 0.85

    def test_raises_plink_error_on_timeout(self, tmp_path, mock_plink_files):
        """Should raise PlinkError when PLINK times out."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.side_effect = subprocess.TimeoutExpired(
                    cmd="plink", timeout=300
                )

                with pytest.raises(PlinkError, match="timed out"):
                    calculate_pairwise_ld(
                        bfile_path=mock_plink_files,
                        snp_list=["rs1", "rs2"],
                        working_dir=str(tmp_path),
                    )


class TestAddSpeciesFlags:
    """Tests for _add_species_flags shared helper."""

    def test_canine_adds_dog_flag(self):
        """Should append --dog for canine species."""
        cmd = ["plink"]
        _add_species_flags(cmd, "canine")
        assert cmd == ["plink", "--dog"]

    def test_dog_alias_adds_the_same_flag(self):
        """The alias the gene layer advertises reaches PLINK too."""
        cmd = ["plink"]
        _add_species_flags(cmd, "dog")
        assert cmd == ["plink", "--dog"]

    def test_feline_adds_chr_set_18(self):
        """Should append --chr-set 18 for feline species."""
        cmd = ["plink"]
        _add_species_flags(cmd, "feline")
        assert cmd == ["plink", "--chr-set", "18"]

    def test_human_adds_no_flags(self):
        """Should not add any flags for human (None species)."""
        cmd = ["plink"]
        _add_species_flags(cmd, None)
        assert cmd == ["plink"]

    def test_species_without_plink_flags_adds_none(self):
        """A species the table does not know needs no chromosome-set flags."""
        cmd = ["plink"]
        _add_species_flags(cmd, "bovine")
        assert cmd == ["plink"]


def _ld_command(species):
    return build_ld_command(
        plink_path="plink",
        bfile_path="data",
        lead_snp="rs1",
        output_path="out",
        species=species,
    )


def _pairwise_command(species):
    return build_pairwise_ld_command(
        plink_path="plink", bfile_path="data", output_path="out", species=species
    )


def _species_flags(cmd):
    """The chromosome-set flags a built command carries, with their values."""
    flags = []
    for index, arg in enumerate(cmd):
        if arg == "--dog":
            flags.append(arg)
        elif arg == "--chr-set":
            flags.extend([arg, cmd[index + 1]])
    return flags


@pytest.mark.parametrize(
    "builder", [_ld_command, _pairwise_command], ids=["ld", "pairwise"]
)
@pytest.mark.parametrize(
    ("species", "flags"),
    [
        ("canine", ["--dog"]),
        ("dog", ["--dog"]),
        ("feline", ["--chr-set", "18"]),
        (None, []),
    ],
)
def test_both_builders_carry_the_same_species_flags(builder, species, flags):
    """Neither builder may drift from the shared species table."""
    assert _species_flags(builder(species)) == flags

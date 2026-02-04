"""Tests for LD calculation module."""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from pylocuszoom.ld import (
    build_ld_command,
    build_pairwise_ld_command,
    calculate_ld,
    calculate_pairwise_ld,
    find_plink,
    parse_ld_output,
    parse_pairwise_ld_output,
)


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

    def test_includes_dog_flag_for_canine_species(self):
        """Command should include --dog for canine species."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            species="canine",
        )
        assert "--dog" in cmd

    def test_includes_chr_set_for_feline_species(self):
        """Command should include --chr-set 18 for feline species."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            species="feline",
        )
        assert "--chr-set" in cmd
        assert "18" in cmd

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

    def test_handles_missing_file(self, tmp_path):
        """Should return empty DataFrame for missing file."""
        result = parse_ld_output(str(tmp_path / "nonexistent.ld"), "rs12345")

        assert len(result) == 0
        assert "SNP" in result.columns
        assert "R2" in result.columns

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

    def test_raises_when_plink_not_found(self, mock_plink_files):
        """Should raise FileNotFoundError when PLINK not found."""
        with patch("pylocuszoom.ld.find_plink", return_value=None):
            with pytest.raises(FileNotFoundError, match="PLINK not found"):
                calculate_ld(
                    bfile_path=mock_plink_files,
                    lead_snp="rs12345",
                )

    def test_returns_empty_dataframe_on_plink_failure(self, tmp_path, mock_plink_files):
        """Should return empty DataFrame when PLINK fails."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1)

                result = calculate_ld(
                    bfile_path=mock_plink_files,
                    lead_snp="rs12345",
                    working_dir=str(tmp_path),
                )

                assert len(result) == 0
                assert "SNP" in result.columns
                assert "R2" in result.columns

    def test_cleans_up_temp_directory(self, mock_plink_files):
        """Should clean up temp directory when working_dir not specified."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1)

                # Get initial temp dir count
                temp_base = tempfile.gettempdir()
                initial_dirs = set(os.listdir(temp_base))

                calculate_ld(
                    bfile_path=mock_plink_files,
                    lead_snp="rs12345",
                    working_dir=None,
                )

                # Check no new dirs remain
                final_dirs = set(os.listdir(temp_base))
                new_dirs = final_dirs - initial_dirs
                snp_scope_dirs = [d for d in new_dirs if d.startswith("snp_scope_ld_")]
                assert len(snp_scope_dirs) == 0

    def test_uses_specified_plink_path(self, tmp_path, mock_plink_files):
        """Should use specified PLINK path instead of auto-detecting."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1)

            calculate_ld(
                bfile_path=mock_plink_files,
                lead_snp="rs12345",
                plink_path="/custom/path/plink",
                working_dir=str(tmp_path),
            )

            # Check the command used the custom path
            call_args = mock_run.call_args
            cmd = call_args[0][0]
            assert cmd[0] == "/custom/path/plink"

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

    def test_includes_dog_flag_for_canine(self):
        """Command should include --dog for canine species."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            species="canine",
        )

        assert "--dog" in cmd

    def test_includes_chr_set_for_feline(self):
        """Command should include --chr-set 18 for feline species."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            species="feline",
        )

        assert "--chr-set" in cmd
        assert "18" in cmd

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

    def test_default_metric_is_r2(self):
        """Command should use --r2 by default."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
        )

        assert "--r2" in cmd

    def test_no_species_flag_for_human(self):
        """Command should not include species flag for human."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            species=None,
        )

        assert "--dog" not in cmd
        assert "--chr-set" not in cmd


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

    def test_returns_empty_for_missing_ld_file(self, tmp_path):
        """Should return empty matrix for missing .ld file."""
        snplist_content = """rs1
rs2"""
        snplist_file = tmp_path / "test.snplist"
        snplist_file.write_text(snplist_content)

        matrix, snp_ids = parse_pairwise_ld_output(
            str(tmp_path / "nonexistent.ld"), str(snplist_file)
        )

        assert matrix.empty
        assert snp_ids == []

    def test_returns_empty_for_missing_snplist_file(self, tmp_path):
        """Should return empty matrix for missing .snplist file."""
        ld_content = """1.0\t0.85
0.85\t1.0"""
        ld_file = tmp_path / "test.ld"
        ld_file.write_text(ld_content)

        matrix, snp_ids = parse_pairwise_ld_output(
            str(ld_file), str(tmp_path / "nonexistent.snplist")
        )

        assert matrix.empty
        assert snp_ids == []

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

    def test_returns_empty_on_plink_failure(self, tmp_path, mock_plink_files):
        """Should return empty DataFrame when PLINK fails."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                matrix, snp_ids = calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2"],
                    working_dir=str(tmp_path),
                )

                assert matrix.empty
                assert snp_ids == []

    def test_writes_snp_list_file(self, tmp_path, mock_plink_files):
        """Should write SNP list to file when snp_list provided."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2", "rs3"],
                    working_dir=str(tmp_path),
                )

                # Check the SNP list file was written
                snp_list_file = tmp_path / "snp_list.txt"
                assert snp_list_file.exists()
                content = snp_list_file.read_text()
                assert "rs1" in content
                assert "rs2" in content
                assert "rs3" in content

    def test_uses_extract_flag_with_snp_list(self, tmp_path, mock_plink_files):
        """Should use --extract flag when snp_list provided."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2"],
                    working_dir=str(tmp_path),
                )

                # Check command included --extract
                call_args = mock_run.call_args
                cmd = call_args[0][0]
                assert "--extract" in cmd

    def test_uses_region_flags(self, tmp_path, mock_plink_files):
        """Should use --chr/--from-bp/--to-bp when region provided."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    chrom=1,
                    start=1000000,
                    end=2000000,
                    working_dir=str(tmp_path),
                )

                # Check command included region flags
                call_args = mock_run.call_args
                cmd = call_args[0][0]
                assert "--chr" in cmd
                assert "--from-bp" in cmd
                assert "--to-bp" in cmd

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

    def test_cleans_up_temp_directory(self, mock_plink_files):
        """Should clean up temp directory when working_dir not specified."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                temp_base = tempfile.gettempdir()
                initial_dirs = set(os.listdir(temp_base))

                calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2"],
                    working_dir=None,
                )

                final_dirs = set(os.listdir(temp_base))
                new_dirs = final_dirs - initial_dirs
                pairwise_dirs = [
                    d for d in new_dirs if d.startswith("snp_scope_pairwise_ld_")
                ]
                assert len(pairwise_dirs) == 0

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

    def test_uses_specified_metric(self, tmp_path, mock_plink_files):
        """Should pass metric parameter to build_pairwise_ld_command."""
        with patch("pylocuszoom.ld.find_plink", return_value="/usr/bin/plink1.9"):
            with patch("subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")

                calculate_pairwise_ld(
                    bfile_path=mock_plink_files,
                    snp_list=["rs1", "rs2"],
                    working_dir=str(tmp_path),
                    metric="dprime",
                )

                call_args = mock_run.call_args
                cmd = call_args[0][0]
                assert "--r" in cmd
                assert "dprime" in cmd

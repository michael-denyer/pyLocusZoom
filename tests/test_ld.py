"""Tests for LD calculation module."""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from pylocuszoom.ld import (
    build_ld_command,
    calculate_ld,
    find_plink,
    parse_ld_output,
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
            mock_which.side_effect = (
                lambda x: "/usr/bin/plink1.9" if x == "plink1.9" else None
            )
            result = find_plink()
            assert result == "/usr/bin/plink1.9"

    def test_falls_back_to_plink(self):
        """Should fall back to plink if plink1.9 not found."""
        with patch("shutil.which") as mock_which:
            mock_which.side_effect = (
                lambda x: "/usr/bin/plink" if x == "plink" else None
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

    def test_includes_dog_flag_for_dog_species(self):
        """Command should include --dog for dog species."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            species="dog",
        )
        assert "--dog" in cmd

    def test_includes_chr_set_for_cat_species(self):
        """Command should include --chr-set 18 for cat species."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            species="cat",
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

"""PLINK file validation, executable discovery and process execution."""

import subprocess
import tempfile
from functools import partial
from unittest.mock import MagicMock, patch

import pytest

from pylocuszoom.exceptions import EmptyLDOutputError, PlinkError, ValidationError
from pylocuszoom.ld import (
    calculate_ld,
    calculate_pairwise_ld,
    find_plink,
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


class TestCalculateLd:
    """Tests for calculate_ld function."""

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


ENTRY_POINTS = [
    pytest.param(partial(calculate_ld, lead_snp="rs12345"), id="calculate_ld"),
    pytest.param(
        partial(calculate_pairwise_ld, snp_list=["rs1", "rs2"]),
        id="calculate_pairwise_ld",
    ),
]
"""Both PLINK entry points share one runner, so they share its failure modes."""


@pytest.mark.parametrize("entry_point", ENTRY_POINTS)
class TestPlinkFailureModes:
    """Each way a PLINK run can fail surfaces the same way from both entry points."""

    def test_raises_when_plink_not_found(self, entry_point, fake_plink):
        bfile, _ = fake_plink

        with patch("pylocuszoom.ld.find_plink", return_value=None):
            with pytest.raises(FileNotFoundError, match="PLINK not found"):
                entry_point(bfile_path=bfile)

    def test_nonzero_exit_raises_plink_error_with_stderr(
        self, entry_point, tmp_path, fake_plink
    ):
        bfile, plink_writes = fake_plink

        with plink_writes(None, returncode=1, stderr="Error: variant not found"):
            with pytest.raises(PlinkError, match="exit code 1") as exc_info:
                entry_point(
                    bfile_path=bfile,
                    plink_path="/usr/bin/plink1.9",
                    working_dir=str(tmp_path),
                )

        assert "variant not found" in str(exc_info.value)

    def test_raises_plink_error_on_timeout(self, entry_point, tmp_path, fake_plink):
        bfile, _ = fake_plink
        timeout = subprocess.TimeoutExpired(cmd="plink", timeout=300)

        with patch("subprocess.run", side_effect=timeout):
            with pytest.raises(PlinkError, match="timed out"):
                entry_point(
                    bfile_path=bfile,
                    plink_path="/usr/bin/plink1.9",
                    working_dir=str(tmp_path),
                )

    def test_cleans_up_temp_directory(
        self, entry_point, tmp_path, monkeypatch, fake_plink
    ):
        """A failed run leaves nothing behind in the directory it created."""
        bfile, plink_writes = fake_plink
        temp_base = tmp_path / "tmpbase"
        temp_base.mkdir()
        monkeypatch.setattr(tempfile, "tempdir", str(temp_base))

        with plink_writes(None, returncode=1, stderr="error"):
            with pytest.raises(PlinkError):
                entry_point(bfile_path=bfile, plink_path="/usr/bin/plink1.9")

        assert list(temp_base.iterdir()) == []


@pytest.mark.parametrize("pairwise", [False, True])
@pytest.mark.parametrize("relative_input", [False, True])
@pytest.mark.parametrize("relative_output", [False, True])
def test_plink_paths_resolve_from_caller_before_changing_directory(
    tmp_path, monkeypatch, pairwise, relative_input, relative_output
):
    import subprocess
    from pathlib import Path

    from pylocuszoom.ld import calculate_pairwise_ld

    monkeypatch.chdir(tmp_path)
    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / f"reference{suffix}").write_bytes(b"")
    executable = tmp_path / "tools" / "plink"
    executable.parent.mkdir()
    executable.write_text("fake executable")

    def child_run(cmd, *, cwd, **kwargs):
        def child_path(arg):
            return Path(cwd) / cmd[cmd.index(arg) + 1]

        assert Path(cwd).is_absolute()
        assert Path(cmd[0]) == executable
        assert child_path("--bfile").with_suffix(".bed").exists()
        out = child_path("--out")
        if pairwise:
            assert child_path("--extract").read_text() == "rs1\nrs2\n"
            Path(f"{out}.ld").write_text("1 0.4\n0.4 1\n")
            Path(f"{out}.snplist").write_text("rs1\nrs2\n")
        else:
            Path(f"{out}.ld").write_text(
                "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n1 100 rs1 1 200 rs2 0.4\n"
            )
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr("subprocess.run", child_run)
    args = dict(
        bfile_path="reference" if relative_input else str(tmp_path / "reference"),
        working_dir="output" if relative_output else str(tmp_path / "output"),
        plink_path="tools/plink" if relative_input else str(executable),
    )
    if pairwise:
        matrix, ids = calculate_pairwise_ld(**args, snp_list=["rs1", "rs2"])
        assert ids == ["rs1", "rs2"]
        assert matrix.loc["rs1", "rs2"] == 0.4
    else:
        result = calculate_ld(**args, lead_snp="rs1")
        assert result.set_index("SNP")["R2"].to_dict() == {"rs1": 1.0, "rs2": 0.4}


def test_bare_executable_name_resolves_relative_path_entry(tmp_path, monkeypatch):
    from pylocuszoom.ld import _resolve_plink

    monkeypatch.chdir(tmp_path)
    executable = tmp_path / "bin" / "plink-custom"
    executable.parent.mkdir()
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o700)
    monkeypatch.setenv("PATH", "bin")
    assert _resolve_plink("plink-custom") == str(executable)


def test_missing_bare_executable_name_raises(tmp_path, monkeypatch):
    from pylocuszoom.ld import _resolve_plink

    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="PLINK not found"):
        _resolve_plink("missing-plink")

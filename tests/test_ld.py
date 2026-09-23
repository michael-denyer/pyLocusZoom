"""PLINK command construction and species flags."""

from unittest.mock import patch

import pytest

from pylocuszoom.exceptions import ValidationError
from pylocuszoom.ld import (
    _add_species_flags,
    build_ld_command,
    build_pairwise_ld_command,
    calculate_pairwise_ld,
)


class TestBuildLdCommand:
    """Tests for build_ld_command function."""

    def test_includes_required_flags(self):
        """Command should include all required PLINK flags."""
        cmd = build_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            lead_snp="rs12345",
            output_path="/path/to/output",
            species="canine",
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
            species="canine",
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
            species="canine",
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
            species="canine",
        )
        assert "--threads" in cmd
        idx = cmd.index("--threads")
        assert cmd[idx + 1] == "4"


class TestBuildPairwiseLdCommand:
    """Tests for build_pairwise_ld_command function."""

    def test_includes_r2_square_flag(self):
        """Command should include --r2 square for pairwise matrix."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            species="canine",
        )

        assert "--r2" in cmd
        assert "square" in cmd

    def test_includes_write_snplist_flag(self):
        """Command should include --write-snplist to track SNP order."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            species="canine",
        )

        assert "--write-snplist" in cmd

    def test_includes_extract_with_snp_list_file(self):
        """Command should include --extract when snp_list_file provided."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            snp_list_file="/path/to/snps.txt",
            species="canine",
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
            species="canine",
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
            species="canine",
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
                species="canine",
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
                    species="canine",
                )
        mock_run.assert_not_called()

    def test_default_metric_is_r2(self):
        """Command should use --r2 by default."""
        cmd = build_pairwise_ld_command(
            plink_path="/usr/bin/plink1.9",
            bfile_path="/path/to/data",
            output_path="/path/to/output",
            species="canine",
        )

        assert "--r2" in cmd


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

    def test_species_without_plink_support_raises(self):
        """Unknown support must not silently use the human chromosome set."""
        cmd = ["plink"]
        with pytest.raises(ValidationError, match="PLINK"):
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


@pytest.mark.parametrize(
    "entry_point",
    [
        lambda: build_ld_command("plink", "data", "rs1", "out"),
        lambda: build_pairwise_ld_command("plink", "data", "out"),
        lambda: calculate_pairwise_ld("data", ["rs1"]),
    ],
    ids=["build_ld_command", "build_pairwise_ld_command", "calculate_pairwise_ld"],
)
def test_species_has_no_default(entry_point):
    """A canine default read human chromosomes 23-26 as dog autosomes silently."""
    with pytest.raises(TypeError, match="species"):
        entry_point()

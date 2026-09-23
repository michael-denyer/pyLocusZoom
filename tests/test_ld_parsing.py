"""Parsing single-lead and pairwise PLINK output files."""

import pytest

from pylocuszoom.exceptions import EmptyLDOutputError, LDUnavailableError, PlinkError
from pylocuszoom.ld import (
    parse_ld_output,
    parse_pairwise_ld_output,
)


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

    def test_rejects_a_variant_id_named_twice(self, tmp_path):
        """A .bim that spells missing ids "." yields output no id can key."""
        ld_file = tmp_path / "l.ld"
        ld_file.write_text(
            "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n"
            "1 1 rs1 1 1 rs1 1\n1 1 rs1 1 5 rs2 0.3\n1 1 rs1 1 6 rs2 0.4\n"
        )

        with pytest.raises(LDUnavailableError, match="rs2"):
            parse_ld_output(str(ld_file), "rs1")


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

    def test_rejects_a_snplist_naming_a_variant_twice(self, tmp_path):
        (tmp_path / "p.ld").write_text("1 0.5\n0.5 1\n")
        (tmp_path / "p.snplist").write_text("rs1\nrs1\n")

        with pytest.raises(LDUnavailableError, match="rs1"):
            parse_pairwise_ld_output(
                str(tmp_path / "p.ld"), str(tmp_path / "p.snplist")
            )
